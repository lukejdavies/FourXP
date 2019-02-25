#' Covert easily-readable common date to Julian date 
#'
#' @description Covert the date in year, month, day and hour in Jullian date
#' 
#' @param year year to convert
#' @param mon month to convert
#' @param mday day to convert
#' @param hour hour to convert
#' @return Julian Date
#' @examples 
#' date2jd(year=2017,mon=12,mday=18,hour=12) 
#' @export
date2jd<-function (year = 1968, mon = 5, mday = 23, hour = 12) 
{
    year[mon %in% 1:2] = year[mon %in% 1:2] - 1
    mon[mon %in% 1:2] = mon[mon %in% 1:2] + 12
    A = floor(year/100)
    B = floor(A/4)
    C = floor(2 - A + B)
    E = floor(365.25 * (year + 4716))
    F = floor(30.6001 * (mon + 1))
    return(C + mday + E + F + hour/24 - 1524.5)
}
