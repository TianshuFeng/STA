% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapper_vis.R
\name{check_color_code}
\alias{check_color_code}
\title{Check the validation of a color code dataframe}
\usage{
check_color_code(color_code)
}
\arguments{
\item{color_code}{The color code dataframe}
}
\description{
\code{check_color_code} check the validation of a color code dataframe
}
\details{
A valid color code dataframe contains two columns: \code{Abbrev}
and \code{Hex.Colors}. The column \code{Abbrev} contains names of
groups that should match the names of groups in the data. The column
\code{Hex.Colors} contains self defined colors of groups where colors should
be defined as Hex color codes. A tool for picking up colors can be found
\href{https://htmlcolorcodes.com/}{here}.
}
\examples{
temp = auto_set_colorcode(c('a','b','c'))
check_color_code(temp)

}
