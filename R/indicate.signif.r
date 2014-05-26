indicate.signif <-
function(x)
{
   symbol <- NULL
   if (x <= 0.1 & x > 0.05) {
      symbol <- "."
      } else if (x <= 0.05 & x > 0.01) {
      symbol <- "*"
      } else if (x <= 0.01 & x > 0.001) {
      symbol <- "**"
      } else if (x <= 0.001) {
      symbol <- "***"
      } else {
      symbol <- " "
   }
   return(symbol)
}
