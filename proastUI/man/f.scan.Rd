% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmdFunctions_data.R
\name{f.scan}
\alias{f.scan}
\title{Load file in proast data format}
\usage{
f.scan(filename = NULL, separator = "")
}
\arguments{
\item{filename}{character string, path to the file with data to be loaded
line 1 with title of the data; line 2 number of variables; line 3 response
type for each variable; line 4 variable names and from line 5 onwards the data}

\item{separator}{character, defines the separator used in the data to be loaded}
}
\value{
List with info, title of the data; nvar, number of variables;
varnames, names of the variables; data, loaded data; dtype, number indicating
the type of response
}
\description{
Load file in proast data format
}

