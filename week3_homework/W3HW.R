if(!require(tidyverse)){
  install.packages("tidyverse")
}

library(tidyverse)
library(pacman)
library("GGally")


#exercise 1.1
missingdata = !is.na(attenu$station)
print(missingdata)
attenu_cleaned = attenu[missingdata,]
head(attenu_cleaned)
dim(attenu_cleaned)

#exercise 1.2
Theoph_2 = Theoph
str(Theoph_2)
med = median(Theoph_2$Dose)
Theoph_2$Dose_Class = ifelse(Theoph_2$Dose >= med, 'high', 'low')
head(Theoph_2)
dim(Theoph_2)

#exercise 1.3
starbucks = read.csv("starbucks.csv")
starbucks
missingdata = !is.na(starbucks)
is_row_empty = rowSums(missingdata)
starbucks_cleaned = starbucks[is_row_empty == 7,]
plot(x = starbucks_cleaned$Calories,y = starbucks_cleaned$Carb,xlab = "Calories"
     ,ylab = "Carbs (grams)")
maxCal = max(starbucks_cleaned$Calories)
drink = starbucks_cleaned[starbucks_cleaned$Calories == maxCal,]
print(drink$Drink)
starbucks_cleaned$is_highest_fat = ifelse(starbucks_cleaned$Fat == 
                                            max(starbucks_cleaned$Fat), 
                                          TRUE, FALSE)
plot(x = starbucks_cleaned$Calories,y = starbucks_cleaned$Carb,
     xlab = "Calories",ylab = "Carbs (grams)", col = 
       factor(starbucks_cleaned$is_highest_fat))

##exercise 1.4
baseball = read.csv("Batting.csv")
print(baseball)
hr3 = baseball[baseball$HR >= 3, ]
dim(hr3)
print (hr3)
plot(baseball$yearID, baseball$HR, xlab = "Year", ylab = "# of HRs")
laa = baseball[baseball$teamID == "LAA",]
plot(laa$yearID, laa$HR, xlab = "Year", ylab = "# of HRs")
atl_pit = baseball[baseball$teamID == "ATL" | baseball$teamID == "PIT",]
plot(atl_pit$yearID, atl_pit$HR, xlab = "Year", ylab = "# of HRs", col = 
       factor(atl_pit$teamID))

  
##exercise 1.5 
easy_plot = function(x, y, color_data){
  mid = median(color_data)
  levels = ifelse(color_data > mid, "low", "high")
  levels = factor(levels)
  plot(x,y, col=factor(levels), pch = 20)
}
easy_plot(starbucks_cleaned$Calories, starbucks_cleaned$Fat, starbucks_cleaned$Fiber)

##exercise 2.1
iris
##Iris dataset gives the measurements in centimeters of variables sepal length and width and petal length and width, respectively, for 50 flowers from each of 3 specifies of iris. The species are Iris setosa, versicolor, and virginica. 
str(iris)
## There are 150 observations for 5 variables. 

##exercise 2.2
print(iris)
## The species variable is categorical while the sepal length, width and petal length and width are continuous. 

##exercise 2.3
hist(iris$Sepal.Length)
hist(iris$Sepal.Width)
## Petal length and width seem to be linearly correlated. 

##exercise 2.4 
avg_sepal_length = mean(iris$Sepal.Width)
iris_copy = iris
comp_sw = ifelse(iris_copy$Sepal.Width > avg_sepal_length, TRUE, FALSE)
iris$csw = comp_sw
boxplot(iris$Sepal.Width ~ iris$csw)

##exercise 2.5
iris = iris[,-6]
pairs(iris)
ggpairs(iris)


#exercise 3.1
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("TCGAbiolinks")
