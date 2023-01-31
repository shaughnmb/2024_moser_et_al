# 2022_moser_et_al

R script and Proteomics Data for:<br/>
**Human iPSC-derived mononuclear phagocytes restore cognition, neural health, and a subpopulation of hippocampal mossy cells in aging mice**<br/>
<br/>
Authors:<br/>
&ensp;&ensp;V. Alexandra Moser<sup>1\*</sup>, Rachel M. Lipman<sup>1</sup>, Shaughn Bell<sup>1</sup>, George Lawless<sup>1</sup>, Luz Jovita Dimas-Harms<sup>1</sup>, Jake Inzalaco<sup>1</sup>, Simion Kreimer<sup>2</sup>, Sarah J. Parker<sup>2</sup>, Helen S. Goodridge<sup>1</sup>, Clive N. Svendsen<sup>1\*\*</sup><br/>
<br/>
Affiliations:<br/>
&ensp;&ensp;<sup>1</sup> Cedars-Sinai Board of Governors Regenerative Medicine Institute, Cedars-Sinai Medical Center, Los Angeles CA, USA<br/>
&ensp;&ensp;<sup>2</sup> Cedars-Sinai Smidt Heart Institute, Department of Cardiology, Cedars-Sinai Medical Center, Los Angeles CA, USA<br/>
&ensp;&ensp;<sup>\*\*</sup> Corresponding author, Email:&ensp;Clive.Svendsen@cshs.org<br/>
<br/>
R version 4.2.2<br/>
R Script Title:&ensp;moser_et_al_2022.R<br/>
R Script Author:&ensp;Shaughn Bell<br/>
R Script Corresponding Email:&ensp;Shaughn.Bell@cshs.org<br/>
<br/>
Notes: <br/>
&ensp;&ensp;A)&ensp;Script makes use of the variables set up under "project information" as<br/>
&ensp;&ensp;&ensp;&ensp;&ensp;well as additional "prefixes" throughout the script for ease of saving<br/>
&ensp;&ensp;&ensp;&ensp;&ensp;files with a similar path and naming structure.  When reloading data <br/>
&ensp;&ensp;&ensp;&ensp;&ensp;(see note "B"), you must either use the full path or reload the prefixes.<br/>
&ensp;&ensp;B)&ensp;Script saves intermediate steps at each major manipulation of the seurat<br/>
&ensp;&ensp;&ensp;&ensp;&ensp;object via the "saveRDS" function.  If needed, these RDS objects can then<br/>
&ensp;&ensp;&ensp;&ensp;&ensp;be reloaded to easily restart at one of these save points without needing<br/> 
&ensp;&ensp;&ensp;&ensp;&ensp;to start from scratch.  However, these are not required for analysis, and<br/>
&ensp;&ensp;&ensp;&ensp;&ensp;they can be skipped to save time and disk space.
