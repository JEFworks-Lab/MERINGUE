#' Make a shiny app to interactively explore data
#'
#' @param pos Position
#' @param gexp Gene expression
#' @param results Dataframe of analysis results
#' @param title Page title
#' @param description Page description
#' @param pal Color ramp palette
#' @param scale Boolean of whether to scale and center expression values
#' @param zlim Color range
#' @param ... Additional plotting parameters
#'
#' @return Shiny app
#'
#' @examples
#' data(mOB)
#' pos <- mOB$pos
#' gexp <- as.matrix(normalizeCounts(mOB$counts, log=FALSE))
#' results <- mOB$results
#' results[, 1:3] <- round(results[, 1:3], digits=3)
#' makeApp(pos, gexp, results, title='mOB', description='Mouse Olfactory Bulb Spatial Transcriptomics Data')
#'
#' @export
#'
makeApp <- function(pos, gexp, results, title, description=NULL,
                    pal=colorRampPalette(c('blue', 'grey', 'red'))(100),
                    scale=TRUE,
                    zlim=c(-1.5,1.5),
                    ...) {

  ## double check
  gexp <- gexp[rownames(results), rownames(pos)]

  if (interactive()) {
    ui <- shiny::fluidPage(
      title = title,
      titlePanel(title),
      shiny::fluidRow(
        column(12, description, hr())
      ),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          DT::dataTableOutput('table')
        ),
        shiny::mainPanel(
          shiny::plotOutput('plot')
        )
      )
    )
    server <- function(input, output) {
      output$table = DT::renderDataTable(results, server = TRUE,
                                      selection = 'single',
                                      options = list(scrollX = TRUE))
      output$plot = renderPlot({
        validate(need(input$table_rows_selected, "Select a Gene (Click on a Data Table Row)"))
        s = input$table_rows_selected
        if(is.null(s)) { s = rownames(gexp)[1] }
        g = rownames(results)[s]
        ggexp = gexp[g, rownames(pos)]
        if(scale) { ggexp = scale(ggexp)[,1] }
        ggexp[ggexp > zlim[2]] <- zlim[2]
        ggexp[ggexp < zlim[1]] <- zlim[1]
        plot(pos, col=map2col(ggexp, pal=pal), pch=16, main=g,
             xlab=NA, ylab=NA, axes=FALSE, ...)
      })
    }
    shiny::shinyApp(ui = ui, server = server)
  }
}


makeMobApp <- function() {
  data(mOB)
  pos <- mOB$pos
  gexp <- as.matrix(normalizeCounts(mOB$counts, log=FALSE))
  results <- mOB$results
  results[, 1:3] <- round(results[, 1:3], digits=3)
  makeApp(pos, gexp, results, title='mOB',
          description='mouse olfactory bulb spatial transcriptomics data',
          cex=3)
}


makeDrosophilaApp <- function() {
  data(drosophila)
  pos <- drosophila$pos
  gexp <- as.matrix(normalizeCounts(drosophila$counts, log=FALSE))
  results <- drosophila$results
  results[, 1:3] <- round(results[, 1:3], digits=3)
  makeApp(pos, gexp, results, title='Drosophila melanogaster embryo',
          description='Drosophila melanogaster embryo aligned in situ hybridization data from BDTN')
}


