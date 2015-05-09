# Simulated p-value
simpval <- 
function(null, obs, alternative)
{
    stopifnot(is.atomic(null))
    stopifnot(is.numeric(null))
    if(length(obs) != 1)
	stop("'obs' must be a vector of length 1!")
    stopifnot(is.numeric(obs))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if(alternative == "two.sided") {
          count <- 2 * min(sum(null <= -abs(obs), na.rm = TRUE),
                sum(null >= abs(obs), na.rm = TRUE))
    } else if(alternative == "less") {
          count <- sum(null <= obs, na.rm = TRUE)
    } else if(alternative == "greater") {
          count <- sum(null >= obs, na.rm=TRUE)
    }
    p <- (count + 1) / (sum(!is.na(null)) + 1)
    return(p)
}
