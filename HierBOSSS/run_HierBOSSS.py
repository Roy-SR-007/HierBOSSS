#############################
# runner file of HierBOSSS
#############################

# required imports

import json
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd
import sympy as sp

from .utils import *

def _regression_metrics(y_true, y_pred):
    """Compute JSON-compatible regression metrics."""
    y_true = np.asarray(y_true, dtype=float).reshape(-1)
    y_pred = np.asarray(y_pred, dtype=float).reshape(-1)

    residual = y_true - y_pred
    ss_res = float(residual @ residual)
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))

    if ss_tot == 0.0:
        r2 = 1.0 if ss_res == 0.0 else None
    else:
        r2 = float(1.0 - ss_res / ss_tot)

    return {
        "RMSE": float(np.sqrt(np.mean(residual ** 2))),
        "MAE": float(np.mean(np.abs(residual))),
        "R2": r2,
    }


def _append_signed_term(
    parts,
    coefficient,
    expression=None,
    significant_digits=6,
):
    """Append one signed term to an expression string."""
    coefficient = float(coefficient)
    magnitude = f"{abs(coefficient):.{significant_digits}g}"

    if expression is None:
        body = magnitude
    else:
        body = f"{magnitude}*({expression})"

    if not parts:
        parts.append(body if coefficient >= 0.0 else f"-{body}")
    else:
        sign = "+" if coefficient >= 0.0 else "-"
        parts.append(f" {sign} {body}")


def _raw_model_string(
    expr_row,
    beta_row,
    K,
    add_intercept,
    significant_digits,
):
    """Combine all raw trees and their posterior coefficients."""
    parts = []

    if add_intercept:
        _append_signed_term(
            parts,
            coefficient=beta_row["intercept"],
            significant_digits=significant_digits,
        )

    for j in range(1, K + 1):
        _append_signed_term(
            parts,
            coefficient=beta_row[f"beta_tree_{j}"],
            expression=str(expr_row[f"tree_{j}"]),
            significant_digits=significant_digits,
        )

    return "".join(parts) if parts else "0"


def _tree_to_sympy(forest, tree_id, node_address=0):
    """Convert a HierBOSSS tree directly into a SymPy expression."""
    node = forest.trees[tree_id][node_address]

    if node is None:
        raise ValueError(
            f"Tree {tree_id + 1} contains an empty node "
            f"at address {node_address}."
        )

    if node.nchild == 0:
        return sp.Symbol(str(node.ft))

    left_address = 2 * node_address + 1
    left = _tree_to_sympy(
        forest=forest,
        tree_id=tree_id,
        node_address=left_address,
    )

    operator = node.op["name"]

    if node.nchild == 1:
        unary_operations = {
            "neg": lambda x: -x,
            "inv": lambda x: 1 / x,
            "sin": sp.sin,
            "cos": sp.cos,
            "exp": sp.exp,
            "sq": lambda x: x ** 2,
            "cu": lambda x: x ** 3,
            "sqrt": lambda x: sp.sqrt(sp.Abs(x)),
        }

        if operator not in unary_operations:
            raise ValueError(
                f"Unsupported unary operator: {operator!r}"
            )

        return unary_operations[operator](left)

    if node.nchild == 2:
        right_address = 2 * node_address + 2
        right = _tree_to_sympy(
            forest=forest,
            tree_id=tree_id,
            node_address=right_address,
        )

        if operator == "add":
            return left + right
        if operator == "mul":
            return left * right

        raise ValueError(
            f"Unsupported binary operator: {operator!r}"
        )

    raise ValueError(f"Unsupported node arity: {node.nchild}")


def _final_expression_and_size(
    detail,
    add_intercept,
    significant_digits,
):
    """
    Construct and simplify the selected expression.

    Model size counts all SymPy nodes, including operators, functions,
    features, coefficient constants, intercepts, and power constants.
    """
    forest = detail["forest"]
    selected_indices = list(detail["selected_indices"])

    beta_full = np.asarray(
        detail["beta_post_mean_full_with_zeros"],
        dtype=float,
    )

    expression = sp.Integer(0)
    tree_offset = int(add_intercept)

    for component_index in selected_indices:
        coefficient = sp.Float(
            str(float(beta_full[component_index])),
            17,
        )

        if add_intercept and component_index == 0:
            component_expression = sp.Integer(1)
        else:
            tree_id = component_index - tree_offset
            component_expression = _tree_to_sympy(
                forest=forest,
                tree_id=tree_id,
            )

        expression += coefficient * component_expression

    # Apply SymPy reductions to the complete selected expression.
    expression = sp.factor(sp.simplify(expression))

    # Round all floating-point constants to the requested precision.
    expression = sp.N(expression, significant_digits)

    model_size = sum(
        1 for _ in sp.preorder_traversal(expression)
    )

    effective_K = sum(
        1
        for component_index in selected_indices
        if not (add_intercept and component_index == 0)
    )

    return (
        str(expression),
        int(model_size),
        int(effective_K),
    )


def _forest_design_matrix(forest, X, add_intercept):
    """Evaluate every symbolic tree and create a design matrix."""
    X = np.asarray(X, dtype=float)

    tree_design = np.column_stack([
        forest._eval_node(0, tree_id, X)
        for tree_id in range(forest.ntrees)
    ])

    if add_intercept:
        design = np.column_stack([
            np.ones(X.shape[0], dtype=float),
            tree_design,
        ])
    else:
        design = tree_design

    if not np.all(np.isfinite(design)):
        raise ValueError(
            "A symbolic tree produced non-finite values."
        )

    return design


def _predict_raw_model(
    forest,
    beta_row,
    X,
    K,
    add_intercept,
):
    """Predict using all trees in the raw forest."""
    design = _forest_design_matrix(
        forest=forest,
        X=X,
        add_intercept=add_intercept,
    )

    coefficients = []

    if add_intercept:
        coefficients.append(float(beta_row["intercept"]))

    coefficients.extend(
        float(beta_row[f"beta_tree_{j}"])
        for j in range(1, K + 1)
    )

    predictions = design @ np.asarray(
        coefficients,
        dtype=float,
    )

    if not np.all(np.isfinite(predictions)):
        raise ValueError(
            "The raw model produced non-finite predictions."
        )

    return predictions


def _predict_final_model(detail, X, add_intercept):
    """Predict using the selected trees and posterior coefficients."""
    design = _forest_design_matrix(
        forest=detail["forest"],
        X=X,
        add_intercept=add_intercept,
    )

    beta_full = np.asarray(
        detail["beta_post_mean_full_with_zeros"],
        dtype=float,
    )

    predictions = design @ beta_full

    if not np.all(np.isfinite(predictions)):
        raise ValueError(
            "The final model produced non-finite predictions."
        )

    return predictions


def _table_metric(metrics, metric_name):
    """Extract a metric for table printing, using NaN for unavailable values."""
    if metrics is None:
        return np.nan

    value = metrics.get(metric_name)
    if value is None:
        return np.nan

    return float(value)


def _supports_ansi_color():
    """Return True when ANSI color codes are appropriate."""
    is_terminal = bool(
        getattr(sys.stdout, "isatty", lambda: False)()
    )
    return (
        is_terminal
        and os.environ.get("NO_COLOR") is None
        and os.environ.get("TERM", "").lower() != "dumb"
    )


def _style_text(text, ansi_code, enabled):
    """Apply optional ANSI styling without adding dependencies."""
    if not enabled:
        return str(text)
    return f"\033[{ansi_code}m{text}\033[0m"


def _result_panel_width():
    """Choose a readable banner width for the current terminal."""
    terminal_width = shutil.get_terminal_size(
        fallback=(100, 24)
    ).columns
    return max(78, min(int(terminal_width), 120))


def _print_results_banner(width, color_enabled):
    """Print the main HierBOSSS results banner."""
    title = f"{HierBOSSS_CONSOLE_LABEL}"
    subtitle = "Summary of Results"

    print()
    print(_style_text("╔" + "═" * width + "╗", "1;36", color_enabled))
    print(
        _style_text(
            "║" + title.center(width) + "║",
            "1;36",
            color_enabled,
        )
    )
    print(
        _style_text(
            "║" + subtitle.center(width) + "║",
            "36",
            color_enabled,
        )
    )
    print(_style_text("╚" + "═" * width + "╝", "1;36", color_enabled))


def _print_section_heading(
    title,
    subtitle,
    width,
    color_enabled,
):
    """Print a titled divider above one result table."""
    label = f" {title} "
    rule_length = max(width - len(label), 2)
    heading = label + "─" * rule_length

    print()
    print(_style_text(heading, "1;36", color_enabled))
    print(f"  {subtitle}")


def _format_runtime(seconds, significant_digits):
    """Format elapsed time in a compact, human-readable form."""
    seconds = float(seconds)

    if seconds < 1.0:
        return (
            f"{seconds * 1000:.{significant_digits}g} ms"
        )

    if seconds < 60.0:
        return f"{seconds:.{significant_digits}g} seconds"

    minutes, remaining_seconds = divmod(seconds, 60.0)
    return (
        f"{int(minutes)} min "
        f"{remaining_seconds:.{significant_digits}g} sec"
    )


def _print_result_tables(output, significant_digits):
    """
    Print an elegant HierBOSSS summary and the raw/final tables.

    Test columns are created only when every returned model contains
    test metrics. When no test set was supplied, no test diagnostic
    columns are included in either table.
    """
    models = output["models"]
    raw_rows = []
    final_rows = []

    test_presence = []

    for model_result in models:
        raw_has_test = (
            model_result["raw_model"].get("test_metrics")
            is not None
        )
        final_has_test = (
            model_result["final_model"].get("test_metrics")
            is not None
        )

        if raw_has_test != final_has_test:
            raise ValueError(
                "Raw and final test-metric availability must match."
            )

        test_presence.append(raw_has_test)

    if test_presence and any(test_presence) and not all(test_presence):
        raise ValueError(
            "Test metrics must be available for either all models "
            "or no models."
        )

    has_test_metrics = bool(
        test_presence and all(test_presence)
    )

    for model_result in models:
        raw = model_result["raw_model"]
        final = model_result["final_model"]

        raw_row = {
            "Rank": model_result["rank"],
            "Raw expression": raw["expression"],
            "Train RMSE": _table_metric(
                raw["train_metrics"], "RMSE"
            ),
            "Train MAE": _table_metric(
                raw["train_metrics"], "MAE"
            ),
            "Train R^2": _table_metric(
                raw["train_metrics"], "R2"
            ),
        }

        if has_test_metrics:
            raw_row.update({
                "Test RMSE": _table_metric(
                    raw["test_metrics"], "RMSE"
                ),
                "Test MAE": _table_metric(
                    raw["test_metrics"], "MAE"
                ),
                "Test R^2": _table_metric(
                    raw["test_metrics"], "R2"
                ),
            })

        # HierBOSSS stores JMP on the log scale.
        raw_row["JMP"] = raw["log_JMP"]
        raw_rows.append(raw_row)

        final_row = {
            "Rank": model_result["rank"],
            "Final expression": final["expression"],
            "Train RMSE": _table_metric(
                final["train_metrics"], "RMSE"
            ),
            "Train MAE": _table_metric(
                final["train_metrics"], "MAE"
            ),
            "Train R^2": _table_metric(
                final["train_metrics"], "R2"
            ),
        }

        if has_test_metrics:
            final_row.update({
                "Test RMSE": _table_metric(
                    final["test_metrics"], "RMSE"
                ),
                "Test MAE": _table_metric(
                    final["test_metrics"], "MAE"
                ),
                "Test R^2": _table_metric(
                    final["test_metrics"], "R2"
                ),
            })

        final_row.update({
            "Effective K": final["effective_K"],
            "Model size": final["model_size"],
        })
        final_rows.append(final_row)

    raw_table = pd.DataFrame(raw_rows)
    final_table = pd.DataFrame(final_rows)

    float_formatter = (
        lambda value: f"{value:.{significant_digits}g}"
    )

    color_enabled = _supports_ansi_color()
    panel_width = _result_panel_width()
    models_returned = int(output["top_r_returned"])
    models_requested = int(output["top_r_requested"])
    criterion = output["selection_criterion"]
    ranking_scale = output["ranking_scale"]

    _print_results_banner(
        width=panel_width,
        color_enabled=color_enabled,
    )

    print()
    print(
        _style_text(
            "Run overview",
            "1;37",
            color_enabled,
        )
    )
    print(
        f"  • Ranked models     : "
        f"{models_returned} returned "
        f"({models_requested} requested)"
    )
    print(f"  • Ranking statistic : {ranking_scale}")
    print(f"  • Reduction method  : {criterion} and SymPy")

    if has_test_metrics:
        print(
            _style_text(
                "  ✓ Evaluation        : "
                "training and held-out test metrics",
                "1;32",
                color_enabled,
            )
        )
    else:
        print(
            _style_text(
                "  ℹ Evaluation        : training metrics only",
                "1;33",
                color_enabled,
            )
        )
        print(
            "    No test set was supplied; test metric "
            "columns and test diagnostics are omitted."
        )

    _print_section_heading(
        title="Raw symbolic models",
        subtitle=(
            ""
        ),
        width=panel_width,
        color_enabled=color_enabled,
    )

    with pd.option_context(
        "display.max_colwidth", None,
        "display.width", None,
        "display.max_columns", None,
    ):
        print(
            raw_table.to_string(
                index=False,
                na_rep="NA",
                float_format=float_formatter,
            )
        )

    _print_section_heading(
        title="Final symbolic models",
        subtitle=(
            ""
        ),
        width=panel_width,
        color_enabled=color_enabled,
    )

    with pd.option_context(
        "display.max_colwidth", None,
        "display.width", None,
        "display.max_columns", None,
    ):
        print(
            final_table.to_string(
                index=False,
                na_rep="NA",
                float_format=float_formatter,
            )
        )

    runtime_text = _format_runtime(
        output["parallel_chains_runtime_seconds"],
        significant_digits=significant_digits,
    )

    print()
    print(
        _style_text(
            "─" * panel_width,
            "1;36",
            color_enabled,
        )
    )
    print(
        _style_text(
            f"✓ {HierBOSSS_CONSOLE_LABEL} analysis complete",
            "1;32",
            color_enabled,
        )
    )
    print(f"  Parallel-chain runtime: {runtime_text}")
    print(
        "  JMP values are reported on the log scale."
    )
    print(
        _style_text(
            "─" * panel_width,
            "1;36",
            color_enabled,
        )
    )

def _print_parallel_launch_banner(
    n_chains,
    maxiter,
    n_trees,
    opset,
):
    """Print an elegant summary before parallel MCMC begins."""
    title = f"Running {HierBOSSS_CONSOLE_LABEL}"

    operator_names = [
        op["name"] if isinstance(op, dict) else str(op)
        for op in opset
    ]

    operator_text = ", ".join(operator_names)

    details = [
        f"Parallel MCMC chains : {n_chains:,}",
        f"MCMC iterations      : {maxiter:,} per chain",
        f"Trees per forest     : {n_trees:,}",
        f"Operator set         : [{operator_text}]",
    ]

    content_width = max(
        56,
        len(title),
        *(len(line) for line in details),
    )

    horizontal_rule = "─" * (content_width + 4)

    print()
    print(f"╭{horizontal_rule}╮")
    print(f"│  {title.center(content_width)}  │")
    print(f"├{horizontal_rule}┤")

    for line in details:
        print(f"│  {line.ljust(content_width)}  │")

    print(f"╰{horizontal_rule}╯")
    print()
    

def run_HierBOSSS(
    X_train,
    y_train,
    K,
    maxdepth,
    seeds,
    prior_params=None,
    add_intercept=True,
    wts_init=None,
    wts_prop=None,
    opset=None,
    ftset=None,
    move_weights=None,
    maxiter=1000,
    burnin=0,
    thin=1,
    n_jobs=None,
    show_progress=True,
    report_every=10,
    r=10,
    X_test=None,
    y_test=None,
    force_intercept=False,
    blr_prior=None,
    prior_variance=10.0,
    rcond=None,
    significant_digits=6,
    print_results=False,
    show_trace_plot=False,
    save_trace_plot=False,
    trace_plot_path=None,
    trace_plot_dpi=300,
):
    """
    Run parallel HierBOSSS chains and return the results as one JSON string.

    Required inputs
    ---------------
    X_train
        Training feature matrix with shape (n_train, p).
    y_train
        Training response vector with shape (n_train,).
    K
        Number of trees in each raw HierBOSSS forest.
    maxdepth
        Maximum tree depth.
    seeds
        Nonempty sequence containing one integer seed per chain.

    Defaults
    --------
    prior_params=None
    add_intercept=True
    wts_init=None
    wts_prop=None
    opset=None
    ftset=None
    move_weights=None
    maxiter=1000
    burnin=0
    thin=1
    n_jobs=None
    show_progress=True
    report_every=1
    r=10
    X_test=None
    y_test=None
    force_intercept=False
    blr_prior=None
    prior_variance=10.0
    rcond=None
    significant_digits=6
    print_results=False

    Returns
    -------
    str
        One JSON string containing the runtime and top-r models.

    Notes
    -----
    BIC is fixed as the subset-selection criterion.
    The runtime measures run_parallel_chains() only.
    """
    X_train = np.asarray(X_train, dtype=float)
    y_train = np.asarray(y_train, dtype=float).reshape(-1)

    if X_train.ndim != 2:
        raise ValueError("X_train must be a 2D array.")
    if X_train.shape[0] != y_train.shape[0]:
        raise ValueError(
            "X_train and y_train must have matching rows."
        )
    if not np.all(np.isfinite(X_train)):
        raise ValueError(
            "X_train must contain only finite values."
        )
    if not np.all(np.isfinite(y_train)):
        raise ValueError(
            "y_train must contain only finite values."
        )

    if not isinstance(r, int) or isinstance(r, bool) or r <= 0:
        raise ValueError("r must be a positive integer.")

    if (
        not isinstance(significant_digits, int)
        or isinstance(significant_digits, bool)
        or significant_digits <= 0
    ):
        raise ValueError(
            "significant_digits must be a positive integer."
        )

    if seeds is None:
        raise ValueError(
            "seeds must be a nonempty sequence of integers."
        )

    seeds = list(seeds)

    if not seeds:
        raise ValueError(
            "seeds must be a nonempty sequence of integers."
        )

    if not all(
        isinstance(seed, (int, np.integer))
        and not isinstance(seed, (bool, np.bool_))
        for seed in seeds
    ):
        raise ValueError("Every chain seed must be an integer.")

    seeds = [int(seed) for seed in seeds]

    if maxiter <= 0:
        raise ValueError("maxiter must be positive.")
    if burnin < 0:
        raise ValueError("burnin must be nonnegative.")
    if thin <= 0:
        raise ValueError("thin must be positive.")
    if burnin + thin > maxiter:
        raise ValueError(
            "At least one MCMC state must be retained; "
            "maxiter must be at least burnin + thin."
        )
    if report_every <= 0:
        raise ValueError("report_every must be positive.")

    test_supplied = X_test is not None or y_test is not None

    if test_supplied and (X_test is None or y_test is None):
        raise ValueError(
            "X_test and y_test must be provided together."
        )

    if test_supplied:
        X_test = np.asarray(X_test, dtype=float)
        y_test = np.asarray(y_test, dtype=float).reshape(-1)

        if X_test.ndim != 2:
            raise ValueError("X_test must be a 2D array.")
        if X_test.shape[1] != X_train.shape[1]:
            raise ValueError(
                "X_test must have the same number of "
                "columns as X_train."
            )
        if X_test.shape[0] != y_test.shape[0]:
            raise ValueError(
                "X_test and y_test must have matching rows."
            )
        if not np.all(np.isfinite(X_test)):
            raise ValueError(
                "X_test must contain only finite values."
            )
        if not np.all(np.isfinite(y_test)):
            raise ValueError(
                "y_test must contain only finite values."
            )

    model = HierBOSSS_MCMC(
        X=X_train,
        y=y_train,
        K=K,
        maxdepth=maxdepth,
        prior_params=prior_params,
        add_intercept=add_intercept,
        wts_init=wts_init,
        wts_prop=wts_prop,
        opset=opset,
        ftset=ftset,
    )
    
    # Print the launch message before timing begins so that
    # parallel_runtime measures parallel-chain computation only.
    if show_progress:
        _print_parallel_launch_banner(
            n_chains=len(seeds),
            maxiter=maxiter,
            n_trees=K,
            opset=model.opset
        )

    # Only parallel-chain execution is timed.
    start_time = time.perf_counter()

    parallel_result = model.run_parallel_chains(
        seeds=seeds,
        move_weights=move_weights,
        maxiter=maxiter,
        burnin=burnin,
        thin=thin,
        n_jobs=n_jobs,
        show_progress=show_progress,
        report_every=report_every,
    )

    parallel_runtime = time.perf_counter() - start_time

    requested_chains = int(
        parallel_result["n_chains_requested"]
    )
    completed_chains = int(
        parallel_result["n_chains_completed"]
    )

    if completed_chains != requested_chains:
        raise RuntimeError(
            f"Only {completed_chains} of {requested_chains} "
            f"chains completed. Errors: "
            f"{parallel_result['errors']}"
        )
        
    trace_plot_saved_path = None

    # Plotting occurs after the parallel-chain timer so
    # displaying or saving the figure is not included in
    # parallel_runtime.
    if show_trace_plot or save_trace_plot:
        (
            trace_figure,
            trace_axis,
            trace_plot_saved_path,
        ) = model.plot_parallel_log_jmp(
            show=show_trace_plot,
            save=save_trace_plot,
            save_path=trace_plot_path,
            dpi=trace_plot_dpi,
        )

    # Ranking and all post-processing occur outside the runtime timer.
    top_expr_df, top_beta_df = (
        model.summarize_parallel_topk_forests(top_k=r)
    )

    if top_expr_df.empty:
        raise RuntimeError(
            "No retained MCMC forests are available."
        )

    reduced_df, reduction_details = reduce_topk_jmp_forests(
        mcmc_obj=model,
        top_expr_df=top_expr_df,
        top_beta_df=top_beta_df,
        top_k=r,
        X=X_train,
        y=y_train,
        criterion="BIC",
        force_intercept=force_intercept,
        blr_prior=blr_prior,
        prior_variance=prior_variance,
        rcond=rcond,
        return_candidate_tables=False,
    )

    ranked_models = []

    for position in range(len(reduced_df)):
        expr_row = top_expr_df.iloc[position]
        beta_row = top_beta_df.iloc[position]
        reduced_row = reduced_df.iloc[position]
        detail = reduction_details[position]
        forest = detail["forest"]

        raw_expression = _raw_model_string(
            expr_row=expr_row,
            beta_row=beta_row,
            K=K,
            add_intercept=add_intercept,
            significant_digits=significant_digits,
        )

        final_expression, model_size, effective_K = (
            _final_expression_and_size(
                detail=detail,
                add_intercept=add_intercept,
                significant_digits=significant_digits,
            )
        )

        raw_train_predictions = _predict_raw_model(
            forest=forest,
            beta_row=beta_row,
            X=X_train,
            K=K,
            add_intercept=add_intercept,
        )
        raw_train_metrics = _regression_metrics(
            y_train,
            raw_train_predictions,
        )

        final_train_predictions = _predict_final_model(
            detail=detail,
            X=X_train,
            add_intercept=add_intercept,
        )
        final_train_metrics = _regression_metrics(
            y_train,
            final_train_predictions,
        )

        raw_test_metrics = None
        final_test_metrics = None

        if test_supplied:
            raw_test_predictions = _predict_raw_model(
                forest=forest,
                beta_row=beta_row,
                X=X_test,
                K=K,
                add_intercept=add_intercept,
            )
            raw_test_metrics = _regression_metrics(
                y_test,
                raw_test_predictions,
            )

            final_test_predictions = _predict_final_model(
                detail=detail,
                X=X_test,
                add_intercept=add_intercept,
            )
            final_test_metrics = _regression_metrics(
                y_test,
                final_test_predictions,
            )

        ranked_models.append({
            "rank": int(beta_row["rank"]),
            "source": {
                "chain_id": int(beta_row["chain_id"]),
                "chain_seed": int(beta_row["chain_seed"]),
                "iteration": int(beta_row["iteration"]),
            },
            "raw_model": {
                "expression": raw_expression,
                "log_JMP": float(beta_row["log_JMP"]),
                "train_metrics": raw_train_metrics,
                "test_metrics": raw_test_metrics,
            },
            "final_model": {
                "expression": final_expression,
                "BIC": float(reduced_row["best_BIC"]),
                "effective_K": effective_K,
                "model_size": model_size,
                "train_metrics": final_train_metrics,
                "test_metrics": final_test_metrics,
            },
        })

    output = {
        "parallel_chains_runtime_seconds": float(
            parallel_runtime
        ),
        "trace_plot_path": (
            str(trace_plot_saved_path)
            if trace_plot_saved_path is not None
            else None
        ),
        "ranking_scale": "log_JMP",
        "selection_criterion": "BIC",
        "significant_digits": int(significant_digits),
        "top_r_requested": int(r),
        "top_r_returned": int(len(ranked_models)),
        "models": ranked_models,
    }

    result_json = json.dumps(
        output,
        indent=2,
        allow_nan=False,
    )

    if print_results:
        _print_result_tables(
            output=output,
            significant_digits=significant_digits,
        )

    return result_json
