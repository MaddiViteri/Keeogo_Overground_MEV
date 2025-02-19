%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 02/18/2025 (MATLAB R2024a)

%This program will process motion capture data collected from Vicon Nexus
%2.13 and Noraxon Ultium Insoles and EMGS using the M4 version of the
%associated software

%Data will be analyzed for Kinematic Parameters (Joint Angles, Center of
%Mass Trajectory, Step Symmetry, Stride and Step Length, Gait Speed, Swing
%and Stance Phase) Spatiotemporal Parameters (Step and Stride time,
%Cadence), and Plantar Pressure Parameters (Pressure Distribution, Center
%of Pressure)

%For this program to work, the data hard drive must be connected to the
%computer. you could change the path to a folder on the local device, but
%for privacy and data storage regulations participant data was stored on a
%separate 128gb sansdisk drive

%Related .m files:
%getDataVicon.m, getDataNoraxon.m, DataNoraxonConcat.m, graphTrialsIndv.m,
%graphTrialsComp

%Requires
%Image Processing Toolbox

%Next Features------------------------------------------------------------%
%change cell arrays to structures
%improve file matching so it always reads the "non" part no matter where it is
%don't graph the NAN values or whatever those weird spikes are
%only use getData functions if no data is loaded in the workspace

%-------------------------------------------------------------------------%

%clear workspace
clc;
close all;
clear;

%GLOBAL VARIABLES

trialnames_TRUE = ["NONKeeogo SS01","NONKeeogo SS02","NONKeeogo FS01","NONKeeogo FS02","Keeogo SS01",...
        "Keeogo SS02","Keeogo FS01","Keeogo FS02","MVC LGas","MVC RGas","MVC LHam","MVC RHam","MVC LQuad","MVC RQuad"];

%insole sizes
insole_sm = struct('length',234,'forefootWidth',89,'usMenSize',"n/a",'usWomenSize',"5-7.5");
insole_m = struct('length',249,'forefootWidth',97,'usMenSize','6-7.5','usWomenSize','8-10.5');
insole_l = struct('length',265,'forefootWidth',97,'usMenSize','8-10.5','usWomenSize','11-13');
insole_lx = struct('length',289,'forefootWidth',99,'usMenSize','11-13','usWomenSize','n/a');

trialInsole = insole_l; %change size here for participants

%import associated tools/libraries/packages
import org.opensim.modeling.*

%import data
%browse for subject file 
path_Data = uigetdir; %harddrive that participant data is stored on

%get vicon data
%[data_Vicon] = getDataVicon (path_Data); fprintf ('\n'); %fixes path (idk why just do it)

%Kinematic Parameters

%Spatiotemporal Paramaters

%Planter Pressure Parameters

%get insole data
[data_Noraxon] = getDataNoraxon (path_Data); fprintf ('\n'); 

%GRAPHING
graphTrialsIndv(trialnames_TRUE,data_Noraxon,trialInsole,path_Data);

graphTrialsComp(trialnames_TRUE,data_Noraxon,trialInsole,path_Data);

    