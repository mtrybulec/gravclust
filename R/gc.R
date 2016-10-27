#' @export
gravClust <- function(x, max.steps = 100, ...) {
    step <- 0
    m <- rep(1, nrow(x))
    names(m) <- 1:nrow(x)
    mm <- outer(m, m, "*")

    x0 <- x
    plot(x, col = "black")

    while (step < max.steps) {
        step <- step + 1

        if (nrow(x) == 1) {
            print(sprintf("ALL MERGED at step %d", step))
            break;
        }

        d <- as.matrix(dist(x, ...))
        diag(d) <- NA
        minD <- min(d, na.rm = TRUE)

        a <- 1 / (d^2)
        a[is.na(a)] <- 0
        amm <- a * mm

        v1 <- outer(x[, 1], x[, 1], "-") * amm
        v2 <- outer(x[, 2], x[, 2], "-") * amm

        v <- cbind(rowSums(v1), rowSums(v2))

        maxV <- max(abs(v), na.rm = TRUE)

        shiftFactor <- 1 / 500

        if (maxV * shiftFactor > minD / 2) {
            minIndexes <- which(d == minD, arr.ind = TRUE)[1, ]

            print(sprintf("merging %s to %s", names(m)[minIndexes[1]], names(m)[minIndexes[2]]))

            x <- x[-minIndexes[1], , drop = FALSE]
            m[minIndexes[2]] <- m[minIndexes[2]] + m[minIndexes[1]]
            m <- m[-minIndexes[1]]
            mm <- outer(m, m, "*")
        } else {
            x <- x - v * shiftFactor

            if (step %% 1 == 0) {
                points(x, col = "red")
            }
        }
    }

    text(x0, labels = 1:nrow(x0))
}

#' @export
runExample <- function() {
    set.seed(102)
    x <- matrix(rnorm(20), ncol = 2)
    gravClust(x, max.steps = 200)
}
