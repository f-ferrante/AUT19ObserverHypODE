# Readme
This code requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. 
Any SDP solver can be used. Here we used SDPT3 freely avaialbe at https://github.com/SQLP/SDPT3

The main code that needs to be exececuted is main.m. The paremeters mu and theta needs to be suitably tuned. The code determines the observer gain and simulates the response of the interconnection of the plant and the observer from the initial condition given in the paper. 
