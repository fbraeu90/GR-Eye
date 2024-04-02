# 3D Homogenized Constrained Mixture Model for Eye Growth and Remodeling - A Model of Staphyloma Formation

The homogenized constrained mixture model is implemented in **Abaqus** through a coupled **User-defined MATerial subroutine ([UMAT](https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/sub/default.htm?startat=ch01s01asb31.html))**. This computational model facilitates the simulation of various Growth & Remodeling (G&R) processes in soft biological tissues within the body. Originally designed for aneurysm growth simulations ([Braeu et al. 2017](https://pubmed.ncbi.nlm.nih.gov/27921189/) and [Braeu et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30413985/)), it has been adapted for cardiac G&R ([Amadeus et al. 2023](https://link.springer.com/article/10.1007/s10237-023-01747-w)).

In this context, we have utilized the homogenized constrained mixture model to simulate ophthalmic G&R processes. Specifically, our focus was on investigating whether G&R theory could provide insights into staphyloma formation resulting from local scleral weakening. We simulated three different G&R scenarios by adjusting the rate of mass turnover (i.e. varying the growth parameter k<sub>&sigma;</sub>) and the type of growth over a duration of 13.7 years (approximately 5000 days).: 
- **Scenario (1)**: A growth parameter of k<sub>&sigma;</sub> = 2∙10<sup>-4</sup> 1/days in combination with volumetric growth in thickness direction
- **Scenario (2)**: A growth parameter of k<sub>&sigma;</sub> = 2∙10<sup>-3</sup> 1/days in combination with volumetric growth in thickness direction
- **Scenario (3)**: A growth parameter of k<sub>&sigma;</sub> = 2∙10<sup>-3</sup> 1/days in combination with mass density growth

To run the different scenarios via the console, use the following Abaqus command (assuming your current directory is the repository): 

<code>abaqus job=GR_Eye_ScenarioX user=umat.f</code>

Note: Replace "X" with the scenario number you wish to execute.

For further details and a derivation of the equations, please refer to our publication: [Braeu et al. 2024](LINK TO ARXIV ARTICLE).
<br/><br/>
<p align="center">
  <img src="https://github.com/fbraeu90/GR-Eye/assets/142971506/67b54640-57ba-48a8-9697-ff64b10d2965" width="700">
</p>
<em>
  Left column: Final shape of the posterior pole of the eye after 13.7 years of G&R. (1), (2), and (3) refer to the simulated G&R scenarios. Right column: Maximum growth stimulus (i.e. the difference between collagen fiber stress and homeostatic stress normalized w.r.t. the homeostatic stress) in the PPS over a time span of 13.7 years. 
  Timepoint 0 represents the onset of G&R processes.
</em>
