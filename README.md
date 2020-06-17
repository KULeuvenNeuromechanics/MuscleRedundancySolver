# solveMuscleRedundancy

## Purpose of the software

The original intent of the provided MATLAB code was to solve the muscle redundancy problem using direct collocation as described in De Groote F, Kinney AL, Rao AV, Fregly BJ. Evaluation of direct collocation optimal control problem formulations for solving the muscle redundancy problem. Annals of Biomedical Engineering (2016). http://link.springer.com/article/10.1007%2Fs10439-016-1591-9. 

From v3.0, there are possibilities to, concurrently with solving the muscle redundancy problem, estimate parameters of the modelled muscle-tendon units by using collected EMG and ultrasound data. Optimal fiber length, tendon slack length and tendon stiffness can be set as free variables within the muscle redundancy problem. Experimentally measured fiber lengths can be tracked (US-tracking), the tracking error is a part of the objective function. Details on this parameter estimation problem can be found in Delabastita et al. 2020 (https://link.springer.com/article/10.1007/s10439-019-02395-x). Collected EMG can either be tracked (EMG-tracking) or imposed exactly (EMG-driven). Details on using EMG data in the parameter estimation can be found in Falisse 2016 (https://ieeexplore.ieee.org/document/7748556). Another important feature is that the user can estimate muscle-tendon parameters over different trials of the same movement or from different movements. This allows to make estimation more reliable. We reckon that for solving the muscle redundancy problem OpenSim Moco (https://www.biorxiv.org/content/10.1101/839381v1) might be a more user-friendly and straightforward alternative. However, our software allows the combination of different trials to estimate muscle-tendon parameters. Another difference is in that we use automatic differentiation, while this is not (yet) enabled in Moco. 


## Structure of the code

The code allows to solve three optimal control problems in the following order: 

1. Generic muscle redundancy problem: Solves the muscle redundancy problem using the provided musculoskeletal model, inverse kinematics and inverse dynamics data. 

2. Muscle parameter estimation: Solves the muscle redundancy problem where the variable space can be extended with user-specified muscle tendon parameters. Depending on availability the user can provide experimental muscle fiber length data to be tracked by simulated fiber lengths and/or EMG data to be tracked by simulated muscle excitations. The problem can be made EMG-driven as well, where the excitations will match the EMG signal up to a specified tolerance.

3. Validation muscle redundancy problem: Solves the muscle redundancy problem using the optimized musculoskeletal model, inverse kinematics and/or inverse dynamics data. The idea of this simulation is to analyse whether the parameters (optimized by the estimation problem) predict musculoskeletal behaviour better than the generic model.

The user is off course free to select any of the three described problems. 

Any questions? Please contact us:
maarten.afschrift@kuleuven.be and tom.vanwouwe@kuleuven.be for questions on the parameter optimization algorithm;  friedl.degroote@kuleuven.be; antoine.falisse@kuleuven.be.

## Installation Instruction

Add the main folder and subfolder to your MATLAB path 

```matlab
addpath(genpath('C/......./solveMuscleRedundancy'))).
```

Several software packages are needed to run the program:

- The OpenSim MATLAB interface is used to generate the inputs to the optimal control problem based on a scaled OpenSim model and the solution of inverse kinematics. To this aim, install OpenSim and set up the OpenSim MATLAB interface (OpenSim: https://simtk.org/frs/?group_id=91, OpenSim API: http://simtk-confluence.stanford.edu:8080/display/OpenSim/Scripting+with+Matlab).
- Casadi is used for nonlinear optimization and algorithmic differentiation (https://web.casadi.org/).

## Main Function

SolveMuscleRedundancy is the main function of this program and is used to solve up to three optimal control problems. 

### Input arguments

#### Required input arguments for SolveMuscleRedundancy

1. **model_path**: directory and filename of the scaled OpenSim model (.osim file). The code should work with any OpenSim model with valid muscle-tendon parameters for which OpenSim's Inverse Dynamics and Muscle Analysis Tools generate reliable results. Note that only the muscle-tendon parameters and not the muscle model specified in the osim-file are used (for details see Muscle model).

2. **time**: 1 x 2 MATLAB array with the initial and final time of the analysis in seconds. Initial and final states influence the optimal controls over a period of about 50 ms at the beginning and end of the time interval over which the optimal control problem is solved. Since in practice the initial and final states are generally unknown, problems should be solved for a time interval containing five additional data points (considering a 100Hz sampling frequency) at the beginning and end of the motion cycle. Those additional data points should not be considered in further analyses. The user should thus not be surprised to observe unrealistically high muscle activation at the beginning of the motion (more details in companion paper).

3. **Misc**: miscellaneous input arguments (matlab structure)

Related to the musculoskeletal model:
   - **Misc.DofNames_Input**  is a cell array specifying for which degrees of freedom you want to solve the muscle redundancy problem. Typically the muscle redundancy problem is solved for one leg at a time (there are no muscles spanning both legs).
   - **Misc.MuscleNames_Input** is a cell array that specifies the muscles to be included when solving the muscle redundancy problem. All muscles that actuate (i.e. have a moment arm with respect to) the degrees of freedom specified in *DofNames_Input* will be selected by default if this array is left empty.

Related to the required input files:
   - **Misc.IKfile**: cell array of filenames of the inverse kinematics solution of different motion trials (.mot file).
   - **Misc.IDfile**: cell array of filenames of the inverse dynamics solution of different motion trials (.sto file).

#### Required input arguments when using EMG data

The following input arguments are required to use EMG data:
   - **Misc.EMGconstr**: boolean to select whether you want to track provided EMG signals.
   - **Misc.EMGfile**: cell array of filenames containing EMG data of different motion trials (.mot file). (can be empty )
   - **Misc.EMGSelection**: cell aray with muscles that are constrained/driven by EMG data.

#### Required input arguments when using Ultrasound data

The following input arguments are required to use ultrasound data:
   - **Misc.UStracking**: boolean to select whether you want to track provided muscle fiber lengths.
   - **Misc.USfile**: cell array of filenames containing fiberlength data of different motion trials (.mot file). The fiber length data is usually measured using ultra-sound (US).
   - **Misc.USSelection**: cell array with muscles used in fiber length tracking.

#### Required input arguments for parameter optimization

The following input arguments are required to optimize parameters:

   - **Misc.Estimate_TendonStiffness**: array with names of muscle from which tendon stiffness will be estimated.

   - **Misc.lb_kT_scaling**: lower bound of the scaling factor that will scale the generic tendon stiffness into the optimized tendon stiffness.

   - **Misc.ub_kT_scaling**: upper bound of the scaling factor that will scale the generic tendon stiffness into the optimized tendon stiffness.	

   - **Misc.Estimate_OptimalFiberLength**: array with names of muscle from which optimal fiber length will be estimated.
   
   - **Misc.lb_lMo_scaling**: lower bound of the scaling factor that will scale the generic optimal fiber length into the optimized optimal fiber length.

   - **Misc.ub_lMo_scaling**: upper bound of the scaling factor that will scale the generic optimal fiber length into the optimized optimal fiber length.

   - **Misc.lb_lTs_scaling**: lower bound of the scaling factor that will scale the generic tendon slack length into the optimized tendon slack length.

   - **Misc.ub_lTs_scaling**: upper bound of the scaling factor that will scale the generic tendon slack length into the optimized tendon slack length.


#### Optional input arguments for SolveMuscleRedundancy

Related to flow control:

   - **Misc.MRSbool**: boolean to select whether you want to solve the generic muscle redundancy problem. This will be used as initial guess in the parameter optimization (default = true).

   - **Misc.ValidationBool**: boolean to select whether you want to solve the validation muscle redundancy problem (default = true).

   - **Misc.PlotBool**: boolean to select whether you want to plot lots of output information of intermediate steps in the script.

Related to parameter optimization:

   - **Misc.Coupled_TendonStiffness**: array with names of muscle from which tendon stiffness will be coupled. This means that the generic tendon stiffnesses of these muscles will be scaled with same variable.

   - **Misc.Coupled_fiber_length**: array with names of muscle from which optimal fiber length will be coupled. This means that the generic fiber lengths of these muscles will be scaled with same variable.  

   - **Misc.Coupled_slack_length**: array with names of muscle from which the optimal slack length will be coupled. This means that the generic tendon slack lengths of these muscles will be scaled with same variable.

Related to weights in objective function:

   - **Misc.wlM**: cost function weighting factor for 'tracking fiber lengths' term.

   - **Misc.wEMG**: cost function weighting factor for 'tracking EMG' term.

   - **Misc.wA**: cost function weighting factor for minimizing effort.

   - **Misc.wTres**: cost function weighting factor for minimizing reserve actuator contribution.

   - **Misc.wVm**: cost function weighting factor for minimizing muscle fiber velocities (term mainly for regularization of the optimization).

Related to lowpass filtering of input data:

   - **Misc.f_cutoff_ID**: cutoff frequency for the butterworth recursive low pass filter applied to the inverse dynamics data (default is 6 Hz).

   - **Misc.f_order_ID**: order of the butterworth recursive low pass filter applied to the inverse dynamics data (default is 6).

   - **Misc.f_cutoff_LMT**: cutoff frequency for the butterworth recursive low pass filter applied to the muscle tendon lengths from the muscle analysis (default is 6 Hz).

   - **Misc.f_order_LMT**: order of the butterworth recursive low pass filter applied to the muscle tendon lengths from the muscle analysis (default is 6).		

   - **Misc.f_cutoff_dM**: cutoff frequency for the butterworth recursive low pass filter applied to the muscle moment arms from the muscle analysis (default is 6 Hz).

   - **Misc.f_order_dM**: order of the butterworth recursive low pass filter applied to the muscle moment arms from the muscle analysis (default is 6).

   - **Misc.f_cutoff_IK**: cutoff frequency for the butterworth recursive low pass filter applied to the inverse kinematics data (default is 6 Hz) when performing the muscle analysis to compute muscle-tendon lengths and moment arms.

   - **Misc.f_order_IK**: order of the butterworth recursive low pass filter applied to the inverse kinematics data (default is 6).

Related to transcription:

   - **Misc.Mesh_Frequency**: number of mesh interval per second (default is 100, but a denser mesh might be required to obtain the desired accuracy especially for faster motions).

Related to nominal parameters model:

   - **Misc.kT**: vector with normalized tendon stiffness for the selected muscles. The order should correspond to *Misc.MuscleNames_Input*. The default value is 35 and a lower value corresponds to a more compliant tendon. The default value will be used when left empty. An example is provided in section *Example Gait10dof18m* to set a different stiffness to the Achilles tendon.
   - **Misc.Set_kT_ByName**: cell array to set the tendon stiffness. The first column is a string wit the name of the muscles, second column is the normalised tendon stiffness.


## Output arguments

We provide all state and control trajectories for the different trials and optimal control problems in one Results structure. Trajectories for different trials and optimal control problems are divided in substructures. All trajectories are interpolated on the mesh points.

1. Time: time vector [s] (dimension: Mesh x 1)

2. MExcitation: muscle excitation [-] (dimension: NMuscles x Mesh)

3. MActivation: muscle activation [-] (dimension: NMuscles x Mesh)

4. RActivation.meshPoints: activation of the reserve actuators [Nm] (dimension: Ndof x Mesh)

5. TForcetilde: normalized tendon force  [-] (dimension: NMuscles x Mesh)

6. TForce: tendon force [N] (dimension: NMuscles x Mesh)

7. Fpe: passive muscle force [N] (dimension: NMuscles x Mesh)

8. lMTinterp: muscle tendon length from muscle analysis[m] (dimension: NMuscles x Mesh)

9. lMtildeopt: normalized muscle fiber length [-] (dimension: NMuscles x Mesh)

10. lm: muscle fiber length [m] (dimension: NMuscles x Mesh)

11. vMtilde: normalized muscle fiber velocity [-] (dimension: NMuscles x Mesh)

12. FMltilde: force-length multiplier [-] (dimension: NMuscles x Mesh)

13. FMvtilde: force-velocity multiplier [-] (dimension: NMuscles x Mesh)

14. Param: scaling factors for the different optimized parameters [-]	

15. Misc: for reference of the settings of the simulation we add the Misc structure to the results. 


## Muscle model

The musculotendon properties are fully described in the supplementary materials of the aforementioned publication. Importantly, only the tendon slack length, optimal muscle fiber length, maximal isometric muscle force, optimal pennation angle and maximal muscle fiber contraction velocity are extracted from the referred OpenSim model. Other properties are defined in the code and can be changed if desired. By default, the activation and deactivation time constants are 15 and 60 ms respectively.


## Examples

### Solve the muscle redundancy problem

In this example (Walking_DeGrooteetal2016), we only solve the muscle redundancy problem without parameter estimation. This means that we try to find the optimal muscle excitations (i.e. controls)  that reconstruct the measured inverse dynamic joint moments with minimal excitations and activations squared. This is similar as in DeGroote 2016 (http://link.springer.com/article/10.1007%2Fs10439-016-1591-9) and was the the main aim of v1 and v2 of this software.

You have to select the opensim model, inverse kinematic solution and inverse dynamic solution.

```matlab
model_path  = fullfile(DataPath,'subject1.osim');
Misc.IKfile = {fullfile(DataPath,'Walking_IK.mot')};
Misc.IDfile = {fullfile(DataPath,'Walking_ID.sto')};
```
You can also select the start and end time of the analysis:

```matlab
time=[0.516 1.95]; % Right stance phase (+50ms beginning and end of time interval, more details see manual and publication)
```
And you have to select the degrees of freedom you want to use in this analysis. Note that the muscle redundancy problem will only be solved for these degrees of freedom. For example if you want to solve only for the right ankle joint:
```matlab
Misc.DofNames_Input={'ankle_angle_r'};    % select the DOFs you want to include in the optimization
```
OR if you want to solve for all degrees of freedom in the ankle, knee and hip joint in this model
```matlab
Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'}; 
```
And select an output directory where you want to save the results.
```matlab
Out_path    = fullfile(MyResultsFolder);                    % folder to store results
```

You can also specify some of the optional input arguments to save the results with a specific outputname, plot the main results automatically in a matlab figure .

```matlab
% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% name output
Misc.OutName = 'Walking_';
```
As an additional optional input argument you can specify that you want to run the muscle redundancy solver. This is true by default
```matlab
Misc.MRSBool = 1;
```

And finally solve the muscle redundancy problem
```matlab
[Results,DatStore] = solveMuscleRedundancy(model_path,time,Out_path,Misc);
```

One important thing to note is that you can select the muscles you want to include in this analysis. For example here, we only select calf muscles and tibialis anterior, all other muscles will be removed from the model. WHen you don't use this input argument or leave this empty, the software will select automatically all muscles that span the selected dofs (Misc.DOfNames_Input).
```matlab
Misc.MuscleNames_Input = {'med_gas_l','lat_gas_l','soleus_l','tib_ant_l'}; % select muscles
```

### Info on parameter estimation

You can estimate the optimal fiber length, tendon slack length and tendon stiffness using EMG-data or ultrasound data. For example if you want to estimate the tendon stiffness of the calf muscles.

```matlab
Misc.Estimate_TendonStifness = {'med_gas_l';'lat_gas_l';'soleus_l'}; % Names of muscles of which tendon stifness is estimated
```

You can also select bounds on maximal deviation of the estimated tendon stiffness from the nominal values (i.e. lw < KOpt/Knominal / ub)

```matlab
Misc.lb_kT_scaling = 0.5; % Lower bound
Misc.ub_kT_scaling = 2; % Upper bound 
```

And finally you can couple the (change in) tendon stiffness of multiple muscles. For example if you assume that the calf muscles share the same tendon

```matlab
Misc.Coupled_TendonStifness = {'med_gas_l';'lat_gas_l';'soleus_l'}; % Couple muscles that should have equal tendon stiffness
```

The same approach is used to estimate optimal fiber length and tendon slack length. Note that that when you want to optimize optimal fiber length, by default you also optimize tendon slack length. Selection of muscles:

```matlab
Estimate_OptimalFiberLength = {'med_gas_l';'soleus_l';'lat_gas_l';'tib_ant_l'}; % Names of muscles of which optimal fiber length is estimated - slack length is estimated for these muscles as well
```
Bounds on the ratio between estimated values and nomimal values
```matlab

Misc.lb_lMo_scaling = 0.1; % Lower bound for scaling optimal fiber length
Misc.ub_lMo_scaling = 2.2; % Upper bound for scaling optimal fiber length
Misc.lb_lTs_scaling = 0.9; % Lower bound for scaling tendon slack length
Misc.ub_lTs_scaling = 1.1; % Upper bound for scaling tendon slack length
```
And coupling of muscles (i.e. equal ratio of change in estimated and nominal muscle properties)

```matlab
Misc.Coupled_fiber_length = {'med_gas_l';'lat_gas_l'}; % Couple muscles that should have equal optimal fiber length
Misc.Coupled_slack_length = {'med_gas_l';'lat_gas_l'}; % Couple muscles that should have equal tendon slack length
```


### EMG-information

In this example (Example_EMGWalking), we use EMG-data to
   1) Constrain the simulated muscle activity based on EMG data in the muscle redundancy problem (Add this)
   2) Estimate muscle-tendon parameters (tendon stiffness, tendon slack length and optimal fiber length) of the calf muscles (EMGDriven_simpleAnkle.m) or multiple lower limb muscles (EMGdriven_LowerLimb.m) using an EMG driven approach. This approach is based on Falisse 2016 (https://ieeexplore.ieee.org/document/7748556).

When using EMG data, you have always have to indicate that you want to use EMG data with a boolean and provide the EMG file (.mot format). This EMG file should contain the processed EMG data (filtered + linear enveloppe). 

```matlab
Misc.EMGconstr  = 1;          % Boolean to select EMG constrained option
Misc.EMGfile = {'C/Path/EMG_gait.mot'}; % path to EMG file
```

You also have to select the muscles that will be driven/constrained by EMG data.
```matlab
Misc.EMGSelection = {'tib_ant_l','lat_gas_l','med_gas_l','soleus_l'};
```

In the preferred case, the names of the muscles in the EMG file and in the model correspond. If this is not the case, you can adapt the names in the .mot file as follows: (note that for example the header of the second collumn in the EMG file will be adapted here to 'bifemlh_r')
```matlab
Misc.EMGheaders = {'Time','bifemlh_r','tib_ant_r','per_long_r','lat_gas_r','bifemsh_r','soleus_r','vas_lat_r','vas_med_r','per_brev_l','tib_ant_l','per_long_l','lat_gas_l','med_gas_l','soleus_l','vas_lat_l','vas_med_l','add_long_l','rect_fem_l','tfl_l','glut_med2_l','bifemsh_l','bifemlh_l','glut_med2_r','rect_fem_r'};
```

The relation between EMG data and simulated muscle excitations can be constrained in two ways. First, you can constrain the optimization variable (s) that scales EMG data to simulated muscle activity. (i.e. S EMG = SimExcitation). In case you normalised your EMG data to MVC measurements, this scale factor should be close to 1.
```matlab
Misc.BoundsScaleEMG = [0.9 1.1];  % maximal value to scale EMG
```

Second, you can select bounds on the deviation between meausred EMG data and simulated muscle excitations (i.e. lower bound <  S EMG - SimExcitation < Upper bound). For example if you want a deviation of 0.1 muscle excitations

```matlab
Misc.EMGbounds  = [-0.1 0.1];     % upper and lower bound for difference between simulated and measured muscle activity
```

Or in the special case of EMG-driven, you can impose this as an equality constraint (S EMG - SimExcitation == 0)

```matlab
Misc.EMGbounds  = 0;     % upper and lower bound for difference between simulated and measured muscle activity
```


As additional settings, you can also drive/constrain multiple muscles based on one signal. For example you can use the signal of the medial gastrocnemius in the excitation of the lateral gastrocnemius.
```matlab
Misc.EMG_MuscleCopies = {'med_gas_l','lat_gas_l'};       %  use gastrocnemius medialis EMG to constrain activity of the lateral gastrocn
```
See "Info on parameter estimation" to combine EMG information with parameter estimation.

### Ultrasound-information

You can use ultrasound information to estimate muscle-tendon properties by tracking fiber lengths in the muscle redundancy problem. We provided an example in the folder "Example_UStracking" based on Delabastita et al. 2020 (https://link.springer.com/article/10.1007/s10439-019-02395-x).

You can input ultrasound data by selecting the right motion file. This file should contain fiber length (in m). Note muscle names in the header of this file should be the same as in the OpenSim model.

```matlab
Misc.USfile = {fullfile(DataPath,'PSF_US_stride2.mot');fullfile(DataPath,'+8_US_stride2.mot');fullfile(DataPath,'+15_US_stride2.mot')};
```
You have to set the boolean for ultrasound tracking and select the muscles you want to use in the tracking.

```matlab
Misc.UStracking  = 1;            % Boolean to select US tracking option
Misc.USSelection = {'med_gas_l'; 'soleus_l'}; % select muscles
```

Finally you can set the weight for tracking ultrasound data in the objective function. 

```matlab
Misc.wlM    = 1;                % weight on tracking fiber length: note that
```
Note that increasing this weight most likely results in "extreme" overfitting. Run the validation tool to investigate this.

```matlab
Misc.ValidationBool = 1;
```

## Release Notes

### version 2.1
From v2.1, CasADi can be used as an alternative to GPOPS-II and ADiGator. CasADi is an open-source tool for nonlinear optimization and algorithmic differentiation (https://web.casadi.org/). Results using CasADi and GPOPS-II are very similar (differences can be attributed to the different direct collocation formulations and scaling). We used CasADi's Opti stack, which is a collection of CasADi helper classes that provides a close correspondence between mathematical NLP notation and computer code (https://web.casadi.org/docs/#document-opti). CasADi is actively maintained and developed, and has an active forum (https://groups.google.com/forum/#!forum/casadi-users).

### version 1.1
From v1.1, an implicit formulation of activation dynamics can be used to solve the muscle redundancy problem. Additionally, by using the activation dynamics model proposed by Raasch et al. (1997), we could introduce a nonlinear change of variables to exactly impose activation dynamics in a continuously differentiable form, omitting the need for a smooth approximation such as described in De Groote et al. (2016). A result of this change of variables is that muscle excitations are not directly accessible during the optimization. Therefore, we replaced muscle excitations by muscle activations in the objective function. This implicit formulation is described in *De Groote F, Pipeleers G, Jonkers I, Demeulenaere B, Patten C, Swevers J, De Schutter J. A physiology based inverse dynamic analysis of human gait: potential and perspectives F. Computer Methods in Biomechanics and Biomedical Engineering (2009).* http://www.tandfonline.com/doi/full/10.1080/10255840902788587. Results from both formulations are very similar (differences can be attributed to the slightly different activation dynamics models and cost functions). However, the formulation with implicit activation dynamics (De Groote et al., (2009)) is computationally faster. This can mainly be explained by the omission of a tanh function in the constraint definition, whose evaluation is computationally expensive when solving the NLP.
