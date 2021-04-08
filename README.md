# Spatial Pattern of AIDS Incidence in China in 2018, and its Association with Education and Healthcare

*For the full contents of analysis report, please visit [our online poster](https://spatialepi563.wordpress.com/2020/11/18/group-16/)*.

# Introduction

Both HIV and AIDS infections have reported an increase in China, 14% rise in new infections in 2018, which are mainly due to the increasing sex transmission in HIV/AIDS. Lack of effective sexually education has been a major barrier of prevent and control HIV. Also, Low access to healthcare recourses leads to riskier health behaviors, which can lead to HIV contraction. Hence, we want to explore:
- Is there province-level spatial heterogeneity or autocorrelation of AIDS incidence rate in China in 2018?
- Is the spatial pattern of AIDS incidence rate associated with the spatial variety of education level and healthcare resources?

# Methods

- Performed Chi-square good-of-fit test to decide if there is province-level spatial heterogeneity of AIDS incidence rate.
- Conducted Moran’s I test and created LISA cluster map to see if there is any spatial dependence/autocorrelation of AIDS incidence rate..
- Designed bivariate maps to compare the patterns of EDUCATION and HEALTHCARE with AIDS incidence rate respectively.
- Fitted two linear a-spatial regression models, which are unconditional mean model (model with only intercept but without covariates) and conditional mean model, to see if education level (percentage of population who has high school degree or above) and healthcare resources (number of healthcare institutions) correlate with the AIDS incidence rate. 

# Main Results

## Spatial Heterogeneity

![alt text](https://github.com/Holin-Chen/Chinese-HIV-Spatial-Analysis/blob/main/plots/SMR%20map.png)

- p-value < 0.01, there is statistically significant evidence that the rate of AIDS varies geographically in China
- South-West and Xinjiang (North-West) had significantly higher incidence rates of AIDS than the country-level average rate (p < 0.05), which means the variation is not due to small sample size
- Sichuan (in orange) had the highest rate in China, which was 17.47 cases per 100,000 person

## Spatial Clustering

![alt text](https://github.com/Holin-Chen/Chinese-HIV-Spatial-Analysis/blob/main/plots/LISA%20Cluster%20map.png)

- High-High Cluster: South-West (Sichuan, Yunnan, Chongqing, Guiyang, Guangxi)
- Low-Low Cluster: North-East (Hebei, Shandong)

## Global Regression

### Estimates of coefficients of the conditional model

| **Variable** | **Estimate** |	**Std Error** | **t value** | **p-value** |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| (Intercept) |	10.71 |	4.45 | 2.41 | 0.0228 |
| Healthcare | -4.79E-06 | 3.68E-05 | -0.13 | 0.8973 |
| Education | -0.2181 |	0.1309 | -1.666 | 0.1069 |

The a-spatial regression model shows that both HEALTHCARE and EDUCATION are not statistically significantly associate with RATE, but we can see the p-value for intercept dramatically increased from 4.2e-6 in crude model to 0.0228 in this conditional model. These changes may indicate that although HEALTHCARE and EDUCATION are insignificantly associate with RATE, they can somewhat explain the variation (heterogeneity) of the province-level incident AIDS rates in China.
