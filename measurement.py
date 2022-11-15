import pandas as pd
import numpy as np
import re

from tqdm import tqdm 

from collections import defaultdict, Counter


def num_there(s):
    '''
    입력 string에 숫자가 포함되었는지 확인
    '''
    return any(i.isdigit() for i in s)


def get_digits(s):
    '''
    입력 string에 포함된 숫자를 반환
    '''
    return int(''.join((char for char in s if char.isdigit())))


def get_not_digit(s):
    '''
    입력 string에서 숫자가 아닌 글자를 반환
    '''
    return ''.join((char for char in s if not char.isdigit())).lower()


def preprocessing_drugname(drugname_str):
    '''
    parsing_drugname_str 수행 전 약물명 표현 전처리 (1. 염 정보 제거)
    '''
    salt_names = [' HCl', 'potassium', '/sodium bicarbonate', ' enteric coated']
    drugname_str_prcsd = re.sub('[(].*[)]', '', drugname_str)
    
    for salt_name in salt_names:
        drugname_str_prcsd = drugname_str_prcsd.replace(salt_name, '')

    return drugname_str_prcsd


def find_dose_from_drugname(drugname_str):
    '''
    drugname_str에서 용량 정보 추출하여 반환

    Arg:
        drugname_str (str): 개별 CDW 처방 기록에서의 '약물명(성분명)' 표현

    Outputs:
        dose_str (str): 용량 정보를 담은 string 표현 전체 (ex. '10mg')
        dose_value (int): 단위를 제외한 용량 정보 (ex. 10)
        dose_unit (str): 용량 단위 정보 (ex. 'mg')
    '''
    
    #아래와 같은 용량 표현을 처리할 수 있어야 함... 수정 필요!
    #For dose_str like '10mg/10mL'
    # dose_values = [get_digits(s) for s in dose_str.split('/')]
    # dose_units = [get_not_digit(s) for s in dose_str.split('/')]
    tokens = drugname_str.split()
    dose_str = [token for token in tokens if num_there(token)][0]

    dose_value = get_digits(dose_str)
    dose_unit = get_not_digit(dose_str)

    return dose_str, dose_value, dose_unit


def normalize_drugform(drugform_str):
    '''
    string 형태의 약물 제형 정보를 정규화하여 반환
    '''
    if 'cap' in drugform_str:
        drugform = 'capsule'
    elif 'tab' in drugform_str:
        drugform = 'tablet'
    elif 'inj' in drugform_str:
        drugform = 'injection'
    else:
        drugform = 'none'

    return drugform


def get_dosinginfo_by_row(presc_record):
    '''
    presc_records의 개별 record 내 처방 정보를 요약

    Args:
        presc_record (Dataframe): presc_records 내 개별 행

    Outputs:
        dosinginfo (dictionary): presc_record 내 약물 처방 정보
        - keys: 'drugname', 'dose_value', 'dose_unit', 'drugform'
         * drugname (str): drugname_str 내 약물명 (소문자)
         * dose_value (int): drugname_str에서 약물의 처방 용량 (단위 제외)
         * dose_unit (str): drugname_str에서 약물의 처방 단위
         * drugform (str): drugname_str에서 약물의 제형 ('capsule', 'tablet', 'injection', 'none' 중 하나)
    '''
    drugname_str= presc_record['약품명(성분명)']
    dosing_per_day = presc_record['1일처방량']
    dosing_days = presc_record['투약일수']

    dosinginfo = {}

    #preprocessing drugname
    drugname_str = preprocessing_drugname(drugname_str)

    #parsing drugname_str and get dose information
    dose_str, dose_value, dose_unit = find_dose_from_drugname(drugname_str)

    #parsing drugname and drugform_str
    drugname = drugname_str.split(dose_str)[0].strip().lower()
    drugform_str = ''.join(drugname_str.split(dose_str)[1:]).strip().lower()

    #normalized drugform infomation
    drugform = normalize_drugform(drugform_str)

    #save dosinginfos_by_row
    dosinginfo['drugname'] = drugname
    dosinginfo['dose_value'] = dose_value
    dosinginfo['dose_unit'] = dose_unit
    dosinginfo['drugform'] = drugform
    dosinginfo['dosing_per_day'] = dosing_per_day
    dosinginfo['dosing_days'] = dosing_days

    return dosinginfo


def classify_drugs_by_group(drug_name, drugs_by_group):
    '''
    해당 약물명이 주어진 drugs_by_group에 따라 어디에 속하는지 파악

    Args:
        drugname (str): drug group을 찾고자 하는 약물명 표현
        drugs_by_group (dictionary): group별 약물 성분 정보
         * keys: drug group명 (ex. '<group>ACE inhibitors' 혹은 'lisinopril')
         * values: drug group에 해당하는 약물명 (ex. 'enalapril; peridopril' 혹은 nan)

    Outputs:
        drug_group_classified (str): drugname이 속하는 drug group 표현 (찾지 못한 경우 'none')
    '''
    drug_group_classified = 'none'

    for group, names in drugs_by_group.items():
        candidate_drug_names = names.lower() if group.startswith('<group>') else group.lower()

        if drug_name in candidate_drug_names.split('; '):
            drug_group_classified = group
            break

    return drug_group_classified


def sum_up_summary_presc_by_drug_groups(summary_presc, target_presc, vict_drug, perp_drug):
    '''
    summary_presc를 target_presc에 따라서 다시 정리하여 반환, vict_drug과 perp_drug의 경우에는 구체적인 투약 정보를 요약하여 반환

    Args:
        summary_presc (dictionary): presc records에서 측정한 환자 개인의 투약 정보
        target_presc (dictionary): : 정리하고자 하는 약물 group 및 group 별 약물 성분 정보
         * keys: drug group명 (ex. '<group>ACE inhibitors' 혹은 'lisinopril')
         * values: drug group에 해당하는 약물명 (ex. 'enalapril; peridopril' 혹은 nan)
        vict_drug (str): DDI의 victim drug 약물명
        perp_drug (str): DDI의 perpetrator drug 약물명 

    Output:
        summary_presc_summed (dictionary): target_presc에 따라서 정리한 summary_presc
    '''
    summary_presc_summed = dict()

    for presc_key, presc_value in summary_presc.items():
        drug_name = presc_key.replace('presc_', '')
        drug_group = classify_drugs_by_group(drug_name, target_presc)

        #when drug is victim drug or perpetrator drug
        if drug_name in [vict_drug, perp_drug]:
            presc_key_rev = 'presc_' + drug_name
            summary_presc_summed[presc_key_rev] = presc_value
            summary_presc_summed[presc_key_rev]['prescribed'] = True

        #else when drug is in target_presc 
        elif drug_group!='none':
            #when drug is in drug group
            if drug_group.startswith('<group>'):
                presc_key_rev = 'presc_' + drug_group
                summary_presc_summed[presc_key_rev] = defaultdict(list)
                summary_presc_summed[presc_key_rev]['prescribed'] = True
                summary_presc_summed[presc_key_rev]['total_presc_days'].append(presc_value['total_presc_days'])
            else:
                presc_key_rev = 'presc_' + drug_name
                summary_presc_summed[presc_key_rev] = presc_value
                summary_presc_summed[presc_key_rev]['prescribed'] = True

    #sum up total_presc_days of drug group
    for presc_key, presc_value in summary_presc_summed.items():
        if '<group>' in presc_key:
            summary_presc_summed[presc_key]['total_presc_days'] = sum(presc_value['total_presc_days'])
            #typing to dictionary
            summary_presc_summed[presc_key] = dict(summary_presc_summed[presc_key])

    return summary_presc_summed


def check_diagnosis(record, target_ICD_code):
    '''
    개별 진단 record에서 특정 ICD 코드의 진단 여부를 확인

    Args:
        record (Dataframe): 진단 여부를 확인할 개별 진단 record (단일행)
        target_ICD_code (str): 진단 여부를 확인할 ICD code 표현
         * 자세한 표현 방식은 trial_setting_0210.xlsx 참고

    Return:
        is_diagnosed (bool): 해당 진단 record에서의 진단 발생 여부
    '''  
    def check_diag_X(diag_ICD_code, ICD_code):
        assert('X' in ICD_code, 'Use check_diag_X only when "X" is contained in ICD_code')
        ICD_code_rev = ICD_code[:-2] if ICD_code[:-1].endswith('.') else ICD_code[:-1]

        return diag_ICD_code.startswith(ICD_code_rev)

    ICD_codes = target_ICD_code.replace(' ', '').split(';')

    is_diagnosed = False

    for ICD_code in ICD_codes:
        if ICD_code.endswith('X') or ICD_code.endswith('x'):
            is_diagnosed = check_diag_X(record['ICD10코드'], ICD_code)
        else:
            is_diagnosed = record['ICD10코드'].startswith(ICD_code)

        if is_diagnosed:
            break
    
    return is_diagnosed


def summarize_basic(patient, all_records, index_date):
    '''
    개별 환자의 basic records에서 covariate/outocome 측정 및 요약

    Args:
        patient (Patient): covariate/outocome 측정을 수행할 대상 환자
        all_records (dictionary of Dataframe): covariate/outocome 측정을 수행할 all_records
        index_date (datetime): covariate/outocome 측정에 활용할 index date
         * incident user design의 경우 perpetrator drug incident use, self-controlled case series design의 경우 object drug incident use로 정의

    Output:
        basic_summary (dictionary): basic records에서 측정할 수 있는 환자의 covariate/outocome 항목별로 저장
    '''
    basic_records = all_records['basic']

    basic_summary = {}

    #basic demographics
    basic_summary['ptnt_id'] = patient.all_records['basic']['환자번호'][0]
    age_in_days = (index_date - patient.all_records['basic']['생년월일'][0]).days
    basic_summary['basic_age'] = round(age_in_days/365, 1)
    basic_summary['basic_sex'] = patient.all_records['basic']['성별'][0]
    basic_summary['basic_height'] = round(np.mean(basic_records['신장(cm)']), 1)
    basic_summary['basic_weight'] = round(np.mean(basic_records['체중(kg)']), 1)
    basic_summary['basic_BMI'] = round(np.mean(basic_records['BMI(kg/m²)']), 1)

    #examination schedule
    basic_summary['exam_days_from_first_exam'] = \
        (index_date - patient.all_records['basic']['최초수진일'][0]).days
    basic_summary['exam_days_to_last_exam'] = \
        (patient.all_records['basic']['최종수진일'][0] - index_date).days
    basic_summary['exam_numb'] = len(basic_records['수진(진료)일'])
    
    try:
        basic_summary['exam_reexam_numb'] = \
            basic_records['진료구분(초/재진/신환)'].value_counts()['재진']
    except KeyError:
        basic_summary['exam_reexam_numb'] = 0

    basic_summary['exam_first_exam_numb'] = \
        basic_summary['exam_numb'] - basic_summary['exam_reexam_numb']
    basic_summary['exam_department'] = dict(Counter(basic_records['수진(퇴원포함)진료과']))

    return basic_summary


def summarize_lab(all_records):
    '''
    개별 환자의 lab records에서 covariate/outocome 측정 및 요약

    Args:
        all_records (dictionary of Dataframe): covariate/outocome 측정을 수행할 all_records

    Output:
        lab_summary (dictionary): lab records에서 측정할 수 있는 환자의 covariate/outocome 항목별로 저장
    '''
    lab_records = all_records['lab']

    lab_summary = {}

    # #HbA1c
    lab_records_with_HbA1c = \
        lab_records[['Hb A1c' in lab_code for lab_code in lab_records['검사명']]]
    lab_rslts_HbA1c = lab_records_with_HbA1c['검사결과-수치값']
    lab_summary['HbA1c_numb'] = len(lab_rslts_HbA1c)
    lab_summary['HbA1c_mean'] = round(np.mean(lab_rslts_HbA1c), 2)
    lab_summary['HbA1c_std'] = round(np.std(lab_rslts_HbA1c), 2)

    # numeric, binary, descriptive lab results and codes
    lab_records_with_numeric = \
        lab_records[[not np.isnan(r) for r in lab_records['검사결과-수치값']]]
    lab_records_with_binary = \
        lab_records[[not np.isnan(r) for r in lab_records['검사결과-음성양성']]]
    lab_records_with_descrp = \
        lab_records[[str(r) != 'nan' for r in lab_records['검사결과']]]
  
    lab_summary['lab_numeric_results'] = lab_records_with_numeric['검사결과-수치값']
    lab_summary['lab_numeric_codes']  = lab_records_with_numeric['검사코드']
    lab_summary['lab_binary_results'] = lab_records_with_binary['검사결과-음성양성']
    lab_summary['lab_binary_codes']  = lab_records_with_binary['검사코드']
    lab_summary['lab_descrp_results'] = lab_records_with_descrp['검사결과']
    lab_summary['lab_descrp_codes']  = lab_records_with_descrp['검사코드']
                        
    return lab_summary


def summarize_presc(all_records):
    '''
    개별 환자의 presc records에서 covariate/outcome 측정 및 요약
     * 개별 환자의 CDW 처방 기록을 받아서 처방 정보를 처방 약물별로 투약 정보를 요약하여 반환

    Arg:
        all_records (dictionary of Dataframe): covariate/outcome 측정을 수행할 all_records

    Return:
        presc_summaries (dictionary): 처방 약물 별 투약 정보를 저장
        - keys: 약물명, values: presc_summary
         presc_summary (dictionary): 개별 처방 약물 별 투약 정보를 저장
         - keys: 'drugname', 'total_adm_dose', 'dose_units', 'total_presc_days', 'drugforms'
          * drugname (list of str): 환자가 투약한 약물명 list
          * total_adm_dose (list of int): 해당 처방 기록 내 drugname별 총 투여 용량
          * dose_units (list of str): 해당 처방 기록 내 drugname별 총 용량 단위
          * total_presc_days (list of int): 해당 처방 기록 내 drugname별 총 복용일
          * drugforms (list of str): 해당 처방 기록 내 drugname별 제형 정보
    '''
    presc_records = all_records['presc']

    presc_summaries = {}

    #get dosinginfo from a single presc_record
    dosinginfos = []
    for _, presc_record in presc_records.iterrows():
        dosinginfo = get_dosinginfo_by_row(presc_record)
        dosinginfos.append(dosinginfo)

    #get unique drugnames prescribed in a patient
    drugnames = set([di['drugname'] for di in dosinginfos])

    for drugname in drugnames:
        presc_summary = {}

        #get dose_units and drugforms of a given drugn
        dose_units = list({di['dose_unit'] 
                      for di in dosinginfos if di['drugname']==drugname})
        drugforms = list({dosinginfo['drugform'] 
                      for di in dosinginfos if di['drugname']==drugname})
        
        #calculate total_adm_dose and total_presc_days
        total_adm_dose= 0
        total_presc_days = 0

        for di in dosinginfos:
            if di['drugname']==drugname:
                total_adm_dose += \
                    di['dose_value'] * di['dosing_per_day'] * di['dosing_days'].days
                total_presc_days += di['dosing_days'].days

        #set presc_summary
        presc_summary['total_adm_dose'] = total_adm_dose
        presc_summary['dose_units'] = dose_units
        presc_summary['total_presc_days'] = total_presc_days
        presc_summary['drugforms'] = drugforms
        #save presc_summary
        presc_summaries['presc_' + drugname] = presc_summary

    return presc_summaries


def summarize_diag(all_records, target_ICDs, index_date, only_diagnosed=True):
    '''
    개별 환자의 diag records에서 covariate/outcome 측정 및 요약
     * 개별 환자의 CDW 진단 기록을 받아서 관심 질병별로 진단 정보를 요약하여 반환

    Args:
        all_records (dictionary of Dataframe): covariate/outcome 측정을 수행할 all_records
        target_ICDs (dictionary): covaraite 측정에 사용할 진단 코드 dictionary
         * keys (str): 진단명 (ex. 'Unspecified diabetes mellitus in childbirth')
         * values (str): key에 해당하는 target_ICD (ex. 'O24.9X')
        index_date (datetime): covariate/outcome 측정에 활용할 index date
         * incident user design의 경우 perpetrator drug incident use, self-controlled case series design의 경우 object drug incident use로 정의
        only_diagnosed (bool): True인 경우 진단을 받은 정보만 저장, False의 경우 전체 정보를 저장

    Output:
        diag_summaries (dictionary): 관심 질병별 진단 정보를 저장
        - keys: 'disease_name', value: diag_summary
         * diag_summary (dictionary): 개별 질병의 진단 정보를 저장
          - keys: 'disease_name', 'diagnosed', 'days_to_diagnosied'
          * disease_name (str): 관심 질병의 질병명
          * diagnosed (bool): all_records 내에서 해당 질병의 진단 여부
          * days_from_first_diag (int): index_date에서 해당 질병의 첫 진단으로부터의 날짜 
            (진단을 받은 적 없는 경우 -1)
    '''
    diag_records = all_records['diag']

    diag_summaries = {}

    for target_disease, target_ICD in target_ICDs.items():
        diag_summary = {}
        #get dates_diag from diag_records
        dates_diag = []
        for _, diag_record in diag_records.iterrows():
            if check_diagnosis(diag_record, target_ICD):
                dates_diag.append(diag_record['진단일자'])

        #set diag_summary
        # diag_summary['dates_diag'] = dates_diag
        # diag_summary['used_target_ICD'] = target_ICD
        diag_summary['diagnosed'] = bool(dates_diag)
        diag_summary['days_from_first_diag'] = (min(dates_diag) - index_date).days \
                                                    if bool(dates_diag) else -1
        #save diag_summary
        if only_diagnosed:
            if diag_summary['diagnosed']:
                diag_summaries['diag_' + target_disease] = diag_summary
        else:
            diag_summaries['diag_' + target_disease] = diag_summary

    return diag_summaries


def summarize_cov(patient, all_records_cov, index_date, 
                  target_ICDs, target_presc_cov, vict_drug, perp_drug, drop=['lab']):
    '''
    개별 환자에서 covariate 측정 및 요약

    Args:
        patient (Patient): covariate 측정을 수행할 대상 환자
        all_records_cov (dictionary of Dataframe): covariate 측정을 수행할 records
        index_date (datetime): covariate 측정에 활용할 index date
         * incident user design의 경우 perpetrator drug incident use, self-controlled case series design의 경우 object drug incident use로 정의
        target_ICDs (dictionary): covaraite 측정에 사용할 진단 코드 dictionary
         * keys (str): 진단명 (ex. 'Unspecified diabetes mellitus in childbirth')
         * values (str): key에 해당하는 target_ICD (ex. 'O24.9X')
        target_presc_cov (dictionary): covariate에서 측정하는 약물 group 및 약물 성분 정보
         * keys: drug group명 (ex. '<group>ACE inhibitors' 혹은 'lisinopril')
         * values: drug group에 해당하는 약물명 (ex. 'enalapril; peridopril' 혹은 nan)
        vict_drug (str): DDI의 victim drug 약물명
        perp_drug (str): DDI의 perpetrator drug 약물명 
        drop (list of str): covaraite 측정 및 요약을 수행하지 않은 records type

    Output:
        cov_summary (dictionary): 환자의 covariate를 항목별로 저장 
         * 세부 key 정의는 summarize_basic, lab, presc, diag 함수 정의 참고
    '''
    cov_summary = {}

    # measure covariates from 'basic' records
    if 'basic' not in drop:
        cov_summary_basic = summarize_basic(patient, all_records_cov, index_date)
        cov_summary.update(cov_summary_basic)

    # measure covariates from 'lab' records
    if 'lab' not in drop:
        cov_summary_lab = summarize_lab(all_records_cov)
        cov_summary.update(cov_summary_lab)

    # measure covariates from 'presc' records
    if 'presc' not in drop:
        cov_summary_presc = summarize_presc(all_records_cov)
        cov_summary_presc_summed = sum_up_summary_presc_by_drug_groups(cov_summary_presc,
                                            target_presc=target_presc_cov, 
                                            vict_drug=vict_drug, 
                                            perp_drug=perp_drug)
        cov_summary.update(cov_summary_presc_summed)

    # measure covariates from 'diag' records
    if 'diag' not in drop:
        cov_summary_diag = summarize_diag(all_records_cov, target_ICDs, index_date)
        cov_summary.update(cov_summary_diag)

    return cov_summary


def summarize_otcm(patient, all_records_otcm, index_date, target_ICDs, only_diagnosed=True):
    '''
    개별 환자에서 outcome 측정 및 요약

    Args:
        patient (Patient): outcome 측정을 수행할 대상 환자
        all_records_otcm (dictionary of Dataframe): outcome 측정을 수행할 records
        index_date (datetime): covariate 측정에 활용할 index date
         * incident user design의 경우 perpetrator drug incident use, self-controlled case series design의 경우 object drug incident use로 정의
        target_ICDs (dictionary): covaraite 측정에 사용할 진단 코드 dictionary
         * keys (str): 진단명 (ex. 'Unspecified diabetes mellitus in childbirth')
         * values (str): key에 해당하는 target_ICD (ex. 'O24.9X')
        only_diagnosed (bool): True인 경우 진단을 받은 정보만 저장, False의 경우 전체 정보를 저장

    Output:
        otcm_summary (dictionary): 환자의 outcome 항목별로 저장 
         * 세부 key 정의는 summarize_diag 함수 정의 참고
    '''
    otcm_summary = {}

    #set patient id from basic records
    otcm_summary['ptnt_id'] = patient.all_records['basic']['환자번호'][0]

    #get diag_summary
    otcm_summary_diag = summarize_diag(all_records=all_records_otcm, 
                                       target_ICDs=target_ICDs,
                                       index_date=index_date, 
                                       only_diagnosed=True)

    otcm_summary.update(otcm_summary_diag)

    return otcm_summary


def merging_summary_into_summaries(summary, summaries):
    '''
    전체 환자군의 summaries dictionary에 개별 환자의 covaraite/outcome summary dictionary 정보를 통합하여 반환

    Args:
        summary (dictionary): 개별 환자의 covariate/outcome summary 정보를 담은 dictionary
        summaries (defaultdict, list): 환자 cohort의 covariate/outcome summary 정보를 담은 defaultdict

    Output:
        summaries_merged (defaultdict, list): 개별 환자 정보를 추가한 summaries
    '''
    def merging_dict_in_dict(key_p, value_p):
        '''
        summaries 내 개별 임상정보가 다시 dictionary로 저장된 경우 key 이름을 업데이트하여 
        '''
        if type(value_p)==dict:
            for key_c, value_c in value_p.items():
                merging_dict_in_dict(key_p + '_' + key_c, value_c)
        else:
            summaries[key_p].append(value_p)

    def max_len(summaries):
        '''
        summaries에 이전까지 기록된 환자 수를 파악하여 반환
        '''
        value_lens = [len(v) for v in summaries.values() if type(v)==list]
        return max(value_lens if value_lens else [0])

    summaries_merged = summaries
    old_keys = list(summaries.keys())
    max_len_old = max_len(summaries)

    #Merging a dictionary of a patient
    for key, value in summary.items():
        merging_dict_in_dict(key, value)

    #pad ['nan'] * max_len at the beginning of the list for a new key
    for key in summaries_merged.keys():
        if key not in old_keys:
            summaries_merged[key] = [float('nan')] * max_len_old + summaries_merged[key]

    #Append nan for missing covarites in a patient
    max_len_new = max_len(summaries_merged)
    for key, value in summaries_merged.items():
        if len(value) < max_len_new:
            summaries_merged[key].append(float('nan'))

    #check whether the lengths of the value lists of summarise_merged are all identical
    for key, value in summaries_merged.items():
        assert len(value)==(max_len_old + 1) , f'Please check the length of {key}'

    return summaries_merged


def get_summaries(patients, target_ICDs_cov, target_ICDs_otcm, 
                  target_presc_cov, vict_drug, perp_drug):
    '''
    환자 cohort에서 covariate/outcome 측정하여 분석을 위한 tidy data 형태로 반환

    Arg:
        patients (list of Patient): covariate/outcome을 측정할 patient cohort 내 전체 Patient
        target_ICDs_cov (list of tuple): covariate 측정을 위한 target_ICDs
        target_ICDs_otcm (list of tuple): outcome 측정을 위한 target_ICDs
        target_presc_cov (dictionary): covariate에서 측정하는 약물 group 및 약물 성분 정보
         * keys: drug group명 (ex. '<group>ACE inhibitors' 혹은 'lisinopril')
         * values: drug group에 해당하는 약물명 (ex. 'enalapril; peridopril' 혹은 nan)
        vict_drug (str): DDI의 victim drug 약물명
        perp_drug (str): DDI의 perpetrator drug 약물명 

    Outputs:
        DF_summaries_cov (Dataframe): tidy data 형태로 변환한 summaries_cov 정보
        DF_summaries_otcm (Dataframe): tidy data 형태로 변환한 summaries_otcm 정보
    '''

    summaries_cov = defaultdict(list)
    summaries_otcm = defaultdict(list)

    for patient in tqdm(patients):
        summary_cov= summarize_cov(patient=patient,
                                   all_records_cov=patient.all_records_otcm_itt,
                                   index_date=patient.index_date_final,
                                   target_ICDs=target_ICDs_cov,
                                   target_presc_cov=target_presc_cov,
                                   vict_drug=vict_drug,
                                   perp_drug=perp_drug)
        summary_otcm = summarize_otcm(patient=patient,
                                      all_records_otcm=patient.all_records_otcm_itt,
                                      index_date=patient.index_date_final,
                                      target_ICDs=target_ICDs_otcm)

        summaries_cov = merging_summary_into_summaries(summary_cov, summaries_cov)
        summaries_otcm = merging_summary_into_summaries(summary_otcm, summaries_otcm)

    #save as DataFrame
    DF_summaries_cov = pd.DataFrame(summaries_cov)
    DF_summaries_otcm = pd.DataFrame(summaries_otcm)
    
    return DF_summaries_cov, DF_summaries_otcm
