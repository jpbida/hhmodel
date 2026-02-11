### Running the App ###
Install the required packages

```
install.packages("shiny")
install.packages("deSolve")
```

Source the model file and start the shiny server

```
source("model.R")
shinyApp(ui, server)
```
