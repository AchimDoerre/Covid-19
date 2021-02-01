#age-specific fatality rates
#dI = 1/100*rep(c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3), 2)


#age- and sex-specific infection fatality rates from
#https://www.medrxiv.org/content/10.1101/2020.08.06.20169722v2.full.pdf
dI = 1/100*c(0.005, 0.010, 0.010, 0.035, 0.045, 0.185, 0.575, 2.345, 5.555,
             0.020, 0.000, 0.005, 0.015, 0.085, 0.355, 1.560, 5.530, 14.000)

#dI = 1/100*rep(c(0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.9, 2.1, 5.2), 2)
#age-specific hospitalisation
ha = 1/100*rep(c(0.1, 0.3, 1.2, 3.2, 4.9, 10.2, 16.6, 24.3, 27.3), 2)
#delta = 1.4*gammaI/100*c(0.1, 0.3, 1.2, 3.2, 4.9, 10.2, 16.6, 24.3, 27.3)
ICUp = 1/100*rep(c(5.0, 5.0, 5.0, 5.0, 6.3, 12.2, 27.4, 43.2, 70.9), 2)



shadecolS = rgb(0.75, 0.75, 0.75)
linecolS = rgb(0.2, 0.2, 0.2)


shadecolE = rgb(0.75, 0.75, 0.75)
linecolE = rgb(0.2, 0.2, 0.2)

shadecolI = rgb(1, 217/255, 179/255)
linecolI = rgb(1, 128/255, 0)

shadecolR = rgb(179/255, 1, 179/255)
linecolR = rgb(0.0, 153/255, 0.0)

linecolH = rgb(0.0, 0.0, 1.0)
shadecolH = rgb(179/255, 209/255, 255/255)


shadecolICU = rgb(217/255, 179/255, 1)
linecolICU = rgb(102/255, 0, 204/255)


shadecolD = rgb(1, 153/255, 153/255)
linecolD = rgb(204/255, 0, 0)
