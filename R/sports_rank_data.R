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
#' \item{player_id}{}
#' \item{player_age}{player age at year.}
#' \item{pa}{plate apperances from previous year.}
#' \item{ab}{at bats from previous year.}
#' \item{k_percentage}{k percentage from previous year.}
#' \item{bb_percentage}{bb percentage from previous year.}
#' \item{batting_avg}{batting average from previous year.}
#' \item{slg_percent}{slg percentage from previous year}
#' \item{exit_velocity_avg}{exit velo from previous year}
#' \item{barrel_batted_rate}{barrel batted rate from previous year}
#' \item{whiff_percent}{whiff percentage from previous year}
#' \item{WAR}{Wins Above Replacement from previous year.}
#' \item{BaseSalary}{Base salary from previous year.}
#' \item{n}{number of years ranked.}
#' }
#' @examples
#'   sports_rank_data
"sports_rank_data"
