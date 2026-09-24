$root=Join-Path (Resolve-Path 'dev/check/phase6-package') 'R'
function Replace-Block($file, $pattern, $replacement) {
  $path=Join-Path $root $file
  $text=[IO.File]::ReadAllText($path)
  $regex=[regex]::new($pattern, [Text.RegularExpressions.RegexOptions]::Singleline)
  if($regex.Matches($text).Count -ne 1){throw "Expected exactly one documentation block in $file"}
  $text=$regex.Replace($text,[Text.RegularExpressions.MatchEvaluator]{param($match) $replacement},1)
  [IO.File]::WriteAllText($path,$text)
}
$nu=@'
#' @param nu One finite positive number describing effective Wishart information.
#'   You must supply it; the fit does not estimate it. It is not the marker
#'   count. See the effective-information section and report sensitivity across
#'   plausible values with \code{\link{terradish_rescale_nu}}.
#' @param gradient
'@
Replace-Block 'wishart_covariance.R' "#' @param nu.*?#' @param gradient" $nu
Replace-Block 'generalized_wishart.R' "#' @param nu.*?#' @param gradient" $nu
$covariance=@'
#' The former full-rank likelihood included the arbitrary site-mean direction.
#' It has been removed because centered genetic covariance contains no
#' information in that direction. The contrast form retains the interpretable
#' comparisons among sites and agrees with covariance-derived distances.
#'
#' \strong{Preparing the response.} Use a coherent covariance construction,
#' such as \code{diagonal = "gower"} for grouped data in landgraph. A
#' \code{"within"} diagonal is on a different scale and is unsuitable for
#' Wishart fitting, even when the matrix happens to be positive definite.
#' Per-locus standardization up-weights rare variants; filter by minor-allele
#' frequency and assess sensitivity. Unequal group sizes create unequal
#' sampling variance that a common diagonal nugget does not automatically fix.
#'
#' @seealso
'@
Replace-Block 'wishart_covariance.R' "#' \\strong\{The role of.*?#' @seealso" $covariance
$distance=@'
#' F\eqn{_{ST}} ratio estimators do not have Wishart degrees of freedom. A
#' passed geometry check is necessary but does not establish this sampling
#' likelihood. Prefer covariance-derived squared distances and document the
#' construction; use \code{\link{mlpe}} for ratio-estimator responses.
#'
#' The function checks
'@
Replace-Block 'generalized_wishart.R' "#' \\strong\{The role of.*?#' The function checks" $distance
