crchtree <- function (formula, data, subset, na.action, weights, offset, cluster, ...) {
    control <- mob_control(...)
    if (control$vcov != "opg") {
        warning("only vcov = \"opg\" supported in lmtree")
        control$vcov <- "opg"
    }
    if (!is.null(control$prune)) {
        if (is.character(control$prune)) {
            control$prune <- tolower(control$prune)
            control$prune <- match.arg(control$prune, c("aic", 
                "bic", "none"))
            control$prune <- switch(control$prune, aic = {
                function(objfun, df, nobs) (nobs[1L] * log(objfun[1L]) + 
                  2 * df[1L]) < (nobs[1L] * log(objfun[2L]) + 
                  2 * df[2L])
            }, bic = {
                function(objfun, df, nobs) (nobs[1L] * log(objfun[1L]) + 
                  log(nobs[2L]) * df[1L]) < (nobs[1L] * log(objfun[2L]) + 
                  log(nobs[2L]) * df[2L])
            }, none = {
                NULL
            })
        }
        if (!is.function(control$prune)) {
            warning("unknown specification of 'prune'")
            control$prune <- NULL
        }
    }
    cl <- match.call(expand.dots = TRUE)
    f <- Formula::Formula(formula)
    if (length(f)[2L] == 1L) {
        attr(f, "rhs") <- c(list(1), attr(f, "rhs"))
        formula[[3L]] <- formula(f)[[3L]]
    }
    else {
        f <- NULL
    }
    m <- match.call(expand.dots = FALSE)
    if (!is.null(f)) 
        m$formula <- formula
    m$fit <- lmfit
    m$control <- control
    if ("..." %in% names(m)) 
        m[["..."]] <- NULL
    m[[1L]] <- as.call(quote(partykit::mob))
    rval <- eval(m, parent.frame())
    rval$info$call <- cl
    class(rval) <- c("crchtree", class(rval))
    return(rval)
}

plot.crchtree <- function (x, terminal_panel = node_bivplot, tp_args = list(), tnex = NULL, drop_terminal = NULL, ...) {
    nreg <- if (is.null(tp_args$which)) 
        x$info$nreg
    else length(tp_args$which)
    if (nreg < 1L & missing(terminal_panel)) {
        plot.constparty(as.constparty(x), tp_args = tp_args, 
            tnex = tnex, drop_terminal = drop_terminal, ...)
    }
    else {
        if (is.null(tnex)) 
            tnex <- if (is.null(terminal_panel)) 
                1L
            else 2L * nreg
        if (is.null(drop_terminal)) 
            drop_terminal <- !is.null(terminal_panel)
        plot.modelparty(x, terminal_panel = terminal_panel, tp_args = tp_args, 
            tnex = tnex, drop_terminal = drop_terminal, ...)
    }
}


