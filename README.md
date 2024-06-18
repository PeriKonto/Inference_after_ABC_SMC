Objective: To generate intensity data using wormsim, by assuming some fixed transmission particles. Afterwards the data generated are imported to the ABC-SMC algorithm to recover the original transmission parameter used

Description:

Under the folder Inference with ABC, there are several subfolders, namely the 01_Code, 02_Output, 03_Source_code, the 04_wormsim-v258Ap25. Under the folder 03_Source_code are the r scrips needed to parse values to the xml

Under the folder 01_Code, the r script Model testing.r is needed to generate intensity specific data, given some previously defined transmission parameters (mbr and k). Also possibilities exist to also include value to the external function but for this example I set it to zero. Them xml file required is the template2MDA_all_6_10_18.xml located under the subfolder 03_Source_code. I allowed 1000 repeated simulations and took the average, on intensity which is to be used on ABC. The model outputs is named as Pruda_Wb_based_data_21018.Rdata. The corresponding files to be used later under Task 2.2 are the syn.obs.csv and the ydata.csv. The former is categorized intensity by level and includes the following headers, "village","year","intens","cases","N","prev" and the latter includes the number of people surveyed by age class.
