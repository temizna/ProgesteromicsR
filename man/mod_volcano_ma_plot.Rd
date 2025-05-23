% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_volcano_ma.R
\name{mod_volcano_ma_plot}
\alias{mod_volcano_ma_plot}
\title{Volcano and MA Plot Module}
\usage{
mod_volcano_ma_plot(input, output, session, res_reactive, filtered_data_rv)
}
\arguments{
\item{input}{Shiny input object}

\item{output}{Shiny output object}

\item{session}{Shiny session object}

\item{res_reactive}{ReactiveVal containing DE results}

\item{filtered_data_rv}{(used for species info)}
}
\value{
None. Renders plots and download handlers
}
\description{
Visualizes differential expression results as a volcano plot and MA plot,
includes labeling of top genes and download options
}
