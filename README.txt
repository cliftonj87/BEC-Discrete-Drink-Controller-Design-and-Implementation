This project was made for BME 411 at IUPUI by Team 9 under Dr. Ken Yoshida, and uses several resources that are cited in more depth elsewhere. This project uses the blood ethanol concentration (BEC) model proposed by Dr. Moore et al. in the paper "Pairing food and drink: A physiological model of blood ethanol levels for a variety of drinking behaviors" (https://www.sciencedirect.com/science/article/pii/S0025556422000013).

This text file should be accessible from a GitHub shared through the final report, and aims to outline the files contained within the GitHub related to the project. Files included will be .m, .slx, and .slxc files meant to be run in MATLab and Simulink.

-=- 12/15/22 -=-




-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

"Copy_of_Proj_v3_wF_TF.slx" and "Copy_of_Proj_v3_wF_TF.slxc" are two files involved in a Simulink model. This model is of the third version of the working plant, alongside branches of the second plant being input into a transfer function block. The transfer function is the estimated transfer function that approximates the second plant of the overall system. Meant to simulate through commands present in "Proj_v3_mat.m", as that contains all the variables and parameters.

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

"Proj_v3_1_mat.m" is a MATLab script that is loaded with the parameters and relationships used in the original model, in addition to the updated parameters used in the transformation from monkey to human. The main use for this script was for the "tfest()" function, which is equipped with a for loop that outputs the number of poles, number of zeroes, and the fit percent into a Excel file. All other functions are in the other script.

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

"Proj_v3_mat.m" is a MATLab script that is loaded with the parameters and relationships used in the original model, in addition to the updated parameters used in the transformation from monkey to humnan. The main use for this script is to run all the simulations for the Simulink models present. The variables header creates all the variables that each model needs, the "Run the simulation" header runs the basic Simulink model that doesn't have food and the Simulink model that includes food, the "Tuning the controller" heading was used to change the PID coefficients and the tuning process to compare to the uncompensated system, the "Brute force" header was efforts made to create a transfer function that would estimate the second plant, the "Plotting the arguments of H(x)" runs the Simulink model equipped with the H(x) arguments outputting to the workspace to plot, the "2 TF estimate models" header ran the two plants created to attempt to linearize the system with the two dominant transfer functions, and the "Running standard responses" header ran a Simulink model equipped with a pulse generator block as the input to characterize the second plant.

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

"Proj_v3_wF.slx" and "Proj_v3_wF.slxc" is the Simulink model that was made after the model without food was verified to be working. This was used in the later parts of the project.

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

"Proj_v3_wF_Copy.slx" and "Proj_v3_wF_Copy.slxc" is the Simulink model that is a copy of "Proj_v3_wF" but has To Workspace blocks that return the values of the arguments of the H(x) functions to be observed.

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

"Proj_v3_wF_Steady.slx" and "Proj_v3_wF_Steady.slxc" is the Simulink model with food that was selectively reduced (using the original paper's equations) into a form that a single transfer function could describe. This was done by observing the H(x) arguments, and was used to try and linearize the system. This should have described the steady state of the system after a disturbance, or the very slow and long-term response. It didn't work, but here's the effort.

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

"Proj_v3_wF_Transient.slx" and "Proj_v3_wF_Transient.slxc" is the Simulink model with food that was selectively reduced (using the original paper's equations) into a form that a single transfer function could describe. This was done by observing the H(x) arguments, and was used to try and linearize the system. This should have described the early time of the system after a disturbance, or the fast or transient part of the system. It didn't work, but here's the effort.

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

"Proj_v3_wF_wC.slx" and "Proj_v3_wF_wC.slxc" is the Simulink model with food and a PI controller that aimed to reduce the BEC of the system using food. This was the final model of the system, and the main point of the final presentation.
