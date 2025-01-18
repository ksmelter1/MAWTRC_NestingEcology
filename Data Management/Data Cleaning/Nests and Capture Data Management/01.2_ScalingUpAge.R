#' Read in captures csv
captures <- read_csv("Data Management/Csvs/raw data/captures.csv") 
captures

#' Filter data to include only hens 
captures <- captures %>%
  dplyr::filter(sex == "F") %>%
  dplyr::select(bandid, sex, age, weight, captyr)
captures

#' (?<=_): Only match is there is an underscore immediately before the number we are trying to extract
#' (\\d{4}): Match exactly 4 digits
#' (?=_) : Only match if the four digits are followed 
nest.data <- nest.data %>%
  dplyr::mutate(nestyr = str_extract(nest_id, "(?<=_)(\\d{4})(?=_)"))
glimpse(nest.data)

nest.data$nestyr <- as.numeric(nest.data$nestyr)
nest.data$captyr <- as.numeric(nest.data$captyr)

#' Create a years since capture column
nest.data <- nest.data %>%
  dplyr::mutate(yrsincecap = nestyr-captyr)
glimpse(nest.data)

#' Assign juvenile as the reference level
nest.data.ready$age <- ifelse(nest.data.ready$age == "A", 1, 
                              ifelse(nest.data.ready$age == "J", 0, NA))


#' Dealing with scaling age ad hoc
#' If the bird is an adult and the years since capture is >1 assign it as an adult
#' If not keep the age as juvenile because turkeys will nest the first year as a juvenile
nest.data.ready$age <- ifelse(nest.data.ready$age == "A" & nest.data.ready$yrsincecap >= 1, "A", nest.data.ready$age)
glimpse(nest.data.ready)

summary(nest.data.ready)
str(nest.data.ready$age)