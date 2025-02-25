# ASSIGNMENT 1: Airfoil Characteristics Using Different Analysis Methods

## Assignment Format
A PDF file report of at max 1500 words (3 pages of text, see how to find the word count in Overleaf or Word). There is otherwise not a page limit, but you should only include content that you deem relevant to explain your findings. You should also fulfill the requirements specified below.  
For each group, only one person should upload the final report.  
The deadline for the report can be found on DTU Learn.

## Report Objectives
The assignment aims to assess the advantages and limitations of using different analysis and design tools for airfoils. The tools to be applied are the **Thin Airfoil Theory**, a **Panel Method**, and **Viscous-Inviscid Interaction** using the open-source software **Xfoil**.  

When you are handed the assignment, you are not expected to be able to solve all the tasks yet. 

## Tasks

### 1. Compare Airfoil Geometries
Plot and compare the geometries of the **NACA 2312 and 2324** airfoils, and the **NACA 4412 and 4424** airfoils in a single plot, including their camber lines.

**Remember to use `axis equal` to ensure that the x and y axes have a 1:1 aspect ratio.**

### 2. Lift Coefficient Comparison
Evaluate and compare the lift coefficient using the following five different methods:
- **Thin Airfoil Theory**
- **Panel Method**
- **Xfoil** (with different boundary layer transition settings):
  - Free transition BL
  - Fixed transition BL  

Generate a **plot of lift coefficient vs. angle of attack** ($-10$ to $15$ degrees) using the above methods for each of the four airfoils (**One plot per airfoil**).

**Make sure that in each plot the lines are readable as they can overlap**

**Discussion:** Discuss what causes the differences between the different methods and check if and why some of those differences are more pronounced for some of the airfoils.

### 3. Pressure Difference Distribution
Evaluate and compare the **pressure difference coefficient** $\Delta C_p$ as a function of $x/c$, where:

$$\Delta C_p = \frac{p_{upper} - p_{lower}}{\frac{1}{2} \rho U^2}$$

for **AoA = 10°**, using:
- **Thin Airfoil Theory**
- **Panel Method**
- **Xfoil** (with the settings below):
  - Free transition BL
  - Fixed transition BL  

Generate **pressure coefficient vs. x/c plots** for each airfoil (**One plot per airfoil**). **Ensure the lines are readable and do not overlap excessively.**

**Discussion:** Discuss what causes the differences between the different methods and check if and why some of those differences are more pronounced for one of the airfoils

### 4. Pressure Coefficient Distribution
Do the same as in Task 3 but for the **dimensionless pressure coefficient** $C_p$:

$$C_p = \frac{p - p_0}{\frac{1}{2} \rho U^2}$$

i.e. plot $C_p$ as function of $x/c$, except for **thin airfoil theory**, where the pressure distribution is not known.

### 5. Drag Coefficient and Polar Plot
For the **Xfoil results from Task 2**, save the drag coefficient and generate **polar plots ($C_l$ vs. $C_d$)** for each airfoil.
**Note: Each plot should contain a curve for each boundary layer transition method.**

**Discussion:** Discuss what causes the differences when using different transition models. Check and explain any differences in the trends between the two airfoils.

From the results fill out the follwing table and **discuss the results in it**:
| Case | Maximum Cl/Cd | AoA at Max Cl/Cd |
|------|--------------|------------------|
| NACA 2312 (Free Transition) | | |
| NACA 2312 (Fixed Transition) | | |
| NACA 2324 (Free Transition) | | |
| NACA 2324 (Fixed Transition) | | |
| NACA 4412 (Free Transition) | | |
| NACA 4412 (Fixed Transition) | | |
| NACA 4424 (Free Transition) | | |
| NACA 4424 (Fixed Transition) | | |

### 6. Effect of Lower Reynolds Number
Analyze how a **significantly lower Reynolds number** (e.g. one order of magnitude lower) would affect the results:
- Would the differences between Xfoil and other methods increase or decrease?
- Would the airfoil performance improve or degrade?
- Which airfoil is more sensitive to Reynolds number changes?

Repeat this analysis for different transition criteria by modifying the **N-parameter** and **fixed transition location**.

**Tip:** You can simply run Xfoil with a different Reynolds number and/or transition criteriaand see what happens.

## Required Figures and Tables
- **Airfoil Shape Comparisons:** Plots comparing the airfoil contour/shape comparing the NACA 2312 and 2324 airfoils, and the NACA 4412 and 4424 airfoils, respectively. Remember to include the chamber lines and to use axis equal to ensure that the x and y axis have a 1:1 aspect ratio. (include camber lines, `axis equal` for 1:1 aspect ratio).
- **Lift Coefficient vs. Angle of Attack:** for the 4 different methods/settings, one plot for each of the four airfoils.
- **Pressure Difference Coefficient vs. x/c:** for the 4 different methods/settings, one plot for each of the four airfoils.
- **Pressure Coefficient vs. x/c:** 3 different methods/settings (all except thin airfoil theory), one plot for each of the four airfoils.
- **$C_l$ vs. $C_d$ Polar Plots:** for all the BL transition cases, one plot for each of the four airfoils.
- **Table Summarizing Maximum $C_l/C_d$ and Corresponding AoA.** (see under task 5)

## Xfoil Settings
When using Xfoil, apply the following settings:
- **Reynolds Number:** $1.5 \times 10^6$ (for all evaluations)
- **Mach Number:** $0.0$ (for all evaluations)

### Boundary Layer Transition Models
When evaluating the different Boundary Layer (BL) models assume free transition (Ncrit = 9) for the lower side of the airfoil, and the following transition criteria for the upper side:
- **Free transition:** $N_{crit} = 9$ (lower side of airfoil)
- **Fixed transition:** $N_{crit} = 9$, with transition fixed at $x_t = 0.1$ (upper side)

---
