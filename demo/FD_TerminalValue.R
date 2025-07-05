# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/FD_TerminalValue.R",echo=TRUE)', or
# Type 'demo(FD_TerminalValue)'
# R6 object
FD <- FiniteDifference$new()
# by name
V <- FD$TerminalValue("Linear")
V <- FD$TerminalValue("Stepped")
V <- FD$TerminalValue("Kinked")
V <- FD$TerminalValue("Mitscherlich")
V <- FD$TerminalValue("Gompertz")
V <- FD$TerminalValue("Logistic")
V <- FD$TerminalValue("Transcendental")
V <- FD$TerminalValue("Yield Index")
# by number
V <- FD$TerminalValue(1)
V <- FD$TerminalValue(2)
V <- FD$TerminalValue(3)
V <- FD$TerminalValue(4)
V <- FD$TerminalValue(5)
V <- FD$TerminalValue(6)
V <- FD$TerminalValue(7)
V <- FD$TerminalValue(8)
# user defined terminal value
mu <- FD$get_oup_params()$mu
x <- FD$get_x_stoch_args()$x
V <- (x-mu)^2
FD$set_x_stoch_args(V=V)
FD$set_V_info(name="Quadratic")
FD$TerminalValue()


