library(shiny)
library(shinydashboard)

server <- function(input,output,session) {
  output$val1<-renderText({input$testinput})
}

#Elements of the UI:
header<-dashboardHeader(title = "Coloring_Test",titleWidth = 280)
sidebar<-dashboardSidebar(width = 280,sidebarMenu(id="sidebar_tabs",
                                                  menuItem("AAA", tabName = "AAA")
))

body<-dashboardBody(title="Main",
                    tabItem(tabName = "Overview",h1("Overview"),
                            fluidPage(
                              box(sliderInput(inputId = "testinput",label="testinput",min=-30,max=20,value=5)),
                              # Code with coloring                       
                              box(title="Output",tags$p(textOutput(outputId="val1", inline=TRUE),style="color:#1E90FF")) 
                              # ,
                              # Code with condition, which does not work
                              # box(title="Output",if(textOutput(outputId="val1")>=0){tags$p("IF",style="color:#1E90FF")}else{tags$p("ELSE",style="color:#ff5733")})  # Condition
                            )                  
                    )
)

ui <- dashboardPage(skin = "black", header, sidebar, body)

shinyApp(ui = ui, server = server)