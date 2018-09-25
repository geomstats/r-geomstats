library("scatterplot3d")

# read csv data, add headers
ur5 <- read.csv(file = "ur5testresult_halfspeed_payload1.6lb_1 (1).csv", header = FALSE)
ur5headers <- scan("ur5testresult_header.csv", sep=',', what = "", quiet = TRUE)
colnames(ur5) = ur5headers


x <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (x)`
y <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (y)`
z <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (z)`

scatterplot3d(gsub("[^0-9\\.]", "", x), y, z)
