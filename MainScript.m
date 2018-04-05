clear all
clc

%Load electricity price curve in half hour resolution
elec_price=(3/1000)*xlsread('ElectricityPrice.xlsx'); 
%Repeats every element 3 times. So from half hour to 10 minute resolution
elec_price=repelem(elec_price, 3); 
%Transposes the column to a row
elec_price= elec_price'; 
%Returns the number of intervals in the scheduling horizon
num_interval = size(elec_price,2); 
%Retrieve data for appliances via a spreadsheet
appliance = xlsread('Appliance.xlsx'); 
%Returns the number of appliances by getting the size of the column of
%appliance data
num_app = size(appliance,1); 
%Initialise the schedule matrix to the size of the scheduling horizon as
%rows and respective appliances as columns
schedule = zeros(num_interval,num_app); 

%For every load, perform the appropriate functions
for app_no = 1:num_app 
    
    %Identify the load type
    %%1 represents NLs
    %%2 represents ILs
    %%3 represents ILs with EM
    %%4 represents TCLs
    load_type = appliance(app_no)
    
    if load_type == 1 %NL
        t = appliance(app_no,2); %time interval duration
        tb = appliance(app_no,3); %allowable start time
        te = appliance(app_no,4); %allowable end time
        pw = appliance(app_no,5); %rated power
        l = appliance(app_no,6); %task duration
        pr = elec_price(tb:te); %price during allowable time
        
        solution = NL_ILP(t, pr, l, pw);
        
    end
    if load_type == 2 %IL
        t = appliance(app_no,2); %time interval duration
        tb = appliance(app_no,3); %allowable start time
        te = appliance(app_no,4); %allowable end time
        pw = appliance(app_no,5); %rated power
        e = appliance(app_no,7); %required energy
        pr = elec_price(tb:te); %price during allowable time
        
        solution = IL_ILP( t, pr, pw, e );
        
    end
    if load_type == 3 %IL with EM
        t = appliance(app_no,2); %time interval duration
        tb = appliance(app_no,3); %allowable start time
        te = appliance(app_no,4); %allowable end time
        pw = appliance(app_no,5); %rated power
        e = appliance(app_no,7); %required energy
        t_off = appliance(app_no,8); %min off time
        theta = appliance(app_no,9); %small constant
        pr = elec_price(tb:te); %price during allowable time
        
        solution = IL_EM_ILP( t, pr, pw, e, t_off, theta);
        solution = solution(tb:te); %to get solution of x variables only
    end
    if load_type == 4 %TCL
        t = appliance(app_no,2); %time interval duration
        tb = appliance(app_no,3); %allowable start time
        te = appliance(app_no,4); %allowable end time
        pw = appliance(app_no,5); %rated power
        c = appliance(app_no,10); %specific heat of water
        m = appliance(app_no,11); %mass of water
        temp_up = appliance(app_no,12); %upper limit of water temperature
        temp_0 = appliance(app_no,13); %initial water temperature
        temp_req = appliance(app_no,14); %required temperature
        %Read array of environment temperature from spreadsheet
        temp_en = xlsread('Appliance.xlsx','Sheet2','A:A');
        %Envirtonment temperature during allowable time
        temp_en = temp_en(tb:te); 
        %Read array of demand temperature
        d = xlsread('Appliance.xlsx','Sheet2','B:B');
        %Demand temperature during allowable time
        d = d(tb:te);
        pr = elec_price(tb:te); %price during allowable time
        
        solution = TCL_ILP(t, pr, pw, c, m, temp_up, temp_0, temp_req, temp_en, d);
        
    end
    
    %Fill in array with power status details over scheduling horizon
    schedule(1:num_interval,app_no) = [zeros(tb-1,1);solution;zeros(num_interval - te,1)];
    
end