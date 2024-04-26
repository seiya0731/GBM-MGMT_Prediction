from radiomics import featureextractor
import pandas as pd
import numpy as np
import SimpleITK as sitk
import os

if __name__ == "__main__":
    workRootDir = 'C:\\Users\\Admin\\Desktop\\test_data\\GBM\\MGMT\\Image\\TCGA-GBM\\'
    dataDir = workRootDir + '1_YuChuLi_TCGA\\'
    maskDir = workRootDir + '3_ROI_nnUNet_TCGA\\'

    extractor = featureextractor.RadiomicsFeatureExtractor()
    extractor.loadParams('radiomicsParas.yaml')
    files = os.listdir(maskDir)
    labels =[1, 2, 3]
    modals = ['t1', 't1c', 't2', 't2f']
################################################################################################################
    #check data-mask matching
    # num=0
    # for file in files:
    #     pName = file[0:-7]
    #     imageName = dataDir + pName + modals[3]+' w.nii'
    #     maskName = maskDir + pName + '.nii.gz'
    #     if os.path.isfile(imageName)==False:
    #         print(pName)
    #         num = num + 1
    #         continue
###############################################################################################################
    #check mask ROI
    # num=0 
    # for file in files:
    #     pName = file[0:-7]
    #     maskName = maskDir + pName + '.nii.gz'
    #     mask = sitk.ReadImage(maskName)
    #     data = sitk.GetArrayFromImage(mask)
    #     test = np.zeros_like(data)
    #     test[data==labels[2]] = 1
    #     count = np.sum(test)
    #     if count<=5:
    #         print(pName)
    #         num=num+1
    #         continue      
######################################################################################################################
    #feature extraction
    print('=feature=extraction=started=========================================')
    for modal in modals:
        for label in labels:
            extractor.settings['label'] = label
            df = pd.DataFrame()
            i=1
            for file in files:
                pName = file[0:-7]
                imageName = dataDir + pName + ' ' + modal + ' w.nii'
                maskName = maskDir + pName + '.nii.gz'
                fileName = pName + '_' +modal
                featureVector = extractor.execute(imageName,maskName)
                df_add = pd.DataFrame.from_dict(featureVector.values()).T
                df_add.columns = featureVector.keys()
                patients_ID=pd.DataFrame({'patients_ID':[pName]})
                df_add=pd.concat([patients_ID,df_add],axis=1)
                df = pd.concat([df,df_add])
                print(str(i)+': '+fileName+'-Label='+str(label))
                i=i+1
            print('+++++++++++++++++++++++++All_subjects_Finished++++++++++++++++++++++++++++++++++++++++++++++')
        ######################################################################################################################    
            outDir = workRootDir +'results_'+modal+'-Label='+str(label)+'.csv'
            df.to_csv(outDir)
        

