###R for biologists
##Irina & Rao, 08/02/2023

##Breakout tasks - MSD COVID dataset
#1) Import the corresponding sheet from Excel file provided (use readxl)
#Group 1 - IFNy; Group 2 - IL-1; Group 3 - IL-2; Group 4 - IL-10

#2) Add a new column "Subtraction" (col6) by subtracting the background unstimulated value (col4) from the stimulated one (col3)

#3) Group the "Subtraction" values by the randomisation group (col5) and find mean, median, standard deviation

#4) Find the maximum and minimum values of col6 in each group (group by col5)

#5) Remove the rows, which have negative values in col6

#6) Save the table as csv (comma-separated) file containing columns 1, 2, 5 and 6 only (without negative values - see #5)
