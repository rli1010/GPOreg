trans.reg=function (formula, data, trans = TRUE) 
{#function for fitting the transformed linear model with the GPO
    call <- match.call()
    formula <- as.formula(formula)
    mf <- model.frame(formula, data)
    Ys <- model.response(mf)
    Zmat <- model.matrix(formula, data)
    n <- nrow(data)
    K.x <- ncol(Zmat) - 1
    pihat <- rowMeans(apply(Ys, 2, rank))/n
    pihat[pihat == 1] <- 0.9999
    pihat[pihat == 0] <- 1e-04
    if (trans) {
        etaihat <- qtruncnorm(pihat, -3, 3)
        g.deriv <- 1/dtruncnorm(etaihat, -3, 3)
    }
    else {
        etaihat <- pihat
        g.deriv <- 1
    }
    fit <- lm(as.formula(paste0("etaihat ~", as.character(formula)[[3]])), 
        data = data)
    coefs <- coef(fit)
    var.naive <- vcov(fit)
    s1i <- Zmat * as.vector((etaihat - Zmat %*% coefs))
    Y.ij1 <- apply(Ys, 2, function(x) matrix(rep(x, n), nrow = n, 
        byrow = T) <= matrix(rep(x, n), nrow = n, byrow = F))
    Y.ij2 <- apply(Ys, 2, function(x) (matrix(rep(x, n), nrow = n, 
        byrow = T) < matrix(rep(x, n), nrow = n, byrow = F)) + 
        1/n)
    Y.ij <- (Y.ij1 + Y.ij2)/2
    ksi.ij <- matrix(rowMeans(Y.ij), nrow = n) - pihat
    s2i.ij <- array(NA, c(n, n, K.x + 1))
    for (p in 1:(K.x + 1)) s2i.ij[, , p] <- Zmat[, p] * g.deriv * 
        ksi.ij
    s2i <- apply(s2i.ij, 3, colMeans)
    negA <- crossprod(Zmat)/n
    negA.inv <- solve(negA)
    Var.B <- crossprod(s1i + s2i)/n^2
    var.ABA <- negA.inv %*% Var.B %*% t(negA.inv)
    z.value <- coefs/sqrt(diag(var.ABA))
    p.value <- 2 * (1 - pnorm(abs(z.value)))
    result <- list(coef = coefs, var = var.ABA, z.value = z.value, 
        p.value = p.value, call = call)
    class(result) <- "GPOreg"
    return(result)
}
