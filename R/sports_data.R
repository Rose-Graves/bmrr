#' Ranking of top MLB players from various new sources with covariates
#'
#' Look up airline names from their carrier codes.
#'
#' @source <https://baseballsavant.mlb.com/statcast_search>
#' @format Data frame with 22 columns and 1127 rows
#' \describe{
#' \item{Player}{player name.}
#' \item{Position}{player position.}
#' \item{Year}{year ranked.}
#' \item{ESPN}{rank by ESPN for year.}
#' \item{CBS}{rank by CBS for year.}
#' \item{MLB}{rank by MLB for year.}
#' \item{BR}{rank by Bleacher Report for year.}
#' \item{Yahoo}{rank by Yahoo for year.}
#' \item{player_id}{unique player id.}
#' \item{player_age}{player age at year.}
#' \item{WAR}{Wins Above Replacement from previous year.}
#' \item{BaseSalary}{Base salary from previous year.}
#' }
#' @examples
#'   sports_data
"sports_data"
