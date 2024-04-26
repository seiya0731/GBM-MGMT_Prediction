
'''图像批量预处理'''
import numpy as np
import os
import SimpleITK as sitk


'''调窗'''
def adjustMethod(data_resampled, w_width, w_center):
    val_min = w_center - (w_width / 2)
    val_max = w_center + (w_width / 2)
    data_adjusted = data_resampled.copy()
    data_adjusted[data_resampled < val_min] = val_min
    data_adjusted[data_resampled > val_max] = val_max
    return data_adjusted

'''处理-1000——700数据'''
def pretreatmentImage(image, max=700, min=-1000):
    image[image < min] = min
    image[image > max] = max
    return image

'''0-1归一化操作'''
def MinMax(image):
    if np.max(image) - np.min(image) != 0:
        image = (image - np.min(image)) / \
            (np.max(image) - np.min(image))  # 0-1归一化
    return image

'''z-scoring归一化操作'''
def zScoring(image):
    image = (image - np.mean(image)) / np.std(image)  # z-scoring归一化
    return image

'''重采样'''
def img_resample(ori_data, new_spacing=[1.0, 1.0, 1.0]):

    original_spacing = ori_data.GetSpacing()
    original_size = ori_data.GetSize()
    transform = sitk.Transform()
    transform.SetIdentity()
    new_shape = [
        int(np.round(original_spacing[0] * original_size[0] / new_spacing[0])),
        int(np.round(original_spacing[1] * original_size[1] / new_spacing[1])),
        int(np.round(original_spacing[2] * original_size[2] / new_spacing[2])),
    ]
    # print("新的形状大小为",new_shape)
    resmaple = sitk.ResampleImageFilter()  # 初始化 #初始化一个图像重采样滤波器resampler
    resmaple.SetInterpolator(sitk.sitkLinear)
    resmaple.SetTransform(transform)
    resmaple.SetDefaultPixelValue(0)
    resmaple.SetOutputSpacing(new_spacing)
    resmaple.SetOutputOrigin(ori_data.GetOrigin())
    resmaple.SetOutputDirection(ori_data.GetDirection())
    resmaple.SetSize(new_shape)
    data = resmaple.Execute(ori_data)
    return data

'''N4偏置场校正'''
def correct_bias(in_file, out_file, image_type=sitk.sitkFloat64):
    """
    Corrects the bias using SimpleITK N4BiasFieldCorrection.
    :param in_file: .nii.gz 文件的输入路径
    :param out_file: .nii.gz 校正后的文件保存路径
    :return: 校正后的nii文件全路径名

    """
    # 使用SimpltITK N4BiasFieldCorrection校正MRI图像的偏置场

    input_image = sitk.ReadImage(in_file, image_type)
    output_image_s = sitk.N4BiasFieldCorrection(input_image, input_image > 0)

    return os.path.abspath(out_file)

'''图像几何参数复制'''
def copy_geometry(image: sitk.Image, ref: sitk.Image):
    image.SetOrigin(ref.GetOrigin())
    image.SetDirection(ref.GetDirection())
    image.SetSpacing(ref.GetSpacing())
    return image


if __name__ == "__main__":
    '''读取数据'''
    old_path = 'C:\\Users\\Admin\\Desktop\\test_data\\GBM\MGMT\\Image\\TCGA-GBM\\0_PeiZhun_TCGA\\'
    new_path = 'C:\\Users\\Admin\\Desktop\\test_data\\GBM\\MGMT\\Image\\TCGA-GBM\\1_YuChuLi_TCGA\\'
    new_spacing = [1.0, 1.0, 5.0]
    image_type = sitk.sitkFloat64
    cases = os.listdir(old_path)
    i=0
    for c in cases:
        img_ori = sitk.ReadImage(os.path.join(old_path, c), image_type)
        img_res = img_resample(img_ori, new_spacing)    #图像重采样
        img_N4 = sitk.N4BiasFieldCorrection(img_res, img_res > 0)   #N4偏置场校正
        data_N4 = sitk.GetArrayFromImage(img_N4)    # 获取图像数组
        data_z = zScoring(data_N4)  # Z-scoring归一化
        img_z = sitk.GetImageFromArray(data_z)  # 图像数组转换成sitk对象
        img_out_itk = copy_geometry(img_z, img_ori) # 图像数组转换成sitk对象
        sitk.WriteImage(img_out_itk, os.path.join(new_path, c))
        if i%4==0:
            print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            print('++++++++++++++++++++++++++'+c[0:-9]+'+++++++++++++++++++++++++++++')
            print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print(str(i)+':  '+c)
        i=i+1
