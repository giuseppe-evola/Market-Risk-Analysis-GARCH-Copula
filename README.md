# Advanced VaR & ES Risk Modeling for Multi-Asset Portfolios: GARCH and Copula Approaches

## Overview
This project implements an advanced market risk model for a trading desk managing a 4-stock portfolio. It covers the entire risk assessment pipeline, from log-return computation to GARCH modeling, VaR & ES estimation, and portfolio-level risk analysis using copulas. The study provides insights into model performance through rigorous backtesting and regulatory capital computation.

## Methodology

### 1. Log-Returns & Descriptive Analysis
- Computed daily log-returns of the four stocks.  
- Analyzed heteroscedasticity & normality using statistical tests.  

### 2. GARCH Modeling
- Estimated and compared different GARCH models 
- Selected the best model per stock based on **BIC** criteria.  

### 3. VaR & Expected Shortfall (ES) Estimation
- Used three different methods to compute in-sample **VaR (95%)** and **ES (97.5%)**:
  - **Parametric (Gaussian) VaR**
  - **Historical Simulation (HS)**
  - **Filtered Historical Simulation (FHS)**  

### 4. Backtesting of Risk Models
- Applied **Kupiec and Christoffersen tests** to assess model performance.  
- Evaluated **VaR exceedances and ES violations**.  

### 5. Portfolio-Level Risk Measures
- Assumed **EUR 1 million exposure per stock**.  
- Conducted dependence analysis via **copula models**.  
- Implemented **rolling-window risk estimation**:
  - **FHS with GARCH filtering**
  - **Monte Carlo Simulation with Copula modeling**
- Computed **daily portfolio VaR (99%)** and **regulatory capital**.
