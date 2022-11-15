import pandas as pd
import numpy as np
import random

from tqdm import tqdm

from collections import Counter
from datetime import timedelta

import settings
from patient import check_presc
from supreme_cdw import SupremeCdw

def get_index_dates_IUD(patient, drug_name, washout):
    '''
    Patient 개체를 받아서 incident user design에 맞추어 index_dates를 반환

    Args:
        patient (Patient): index_date를 확인할 Patient 개체
        drug_name (str): incident use를 확인할 약물명
        washout (int): incident use 이전에 확보해야 하는 washout period (단위: days)

    Output:
        index_dates (list of datetime): incident use가 발생한 index_date를 list 형태로 반환
    '''
    index_dates = []
    washout = timedelta(days=washout)

    #set last_adm_date as date of frist examination
    try:
        last_adm_date = patient.all_records['basic']['최초수진일'][0]
    except KeyError:
        print(patient.id)
        print(patient.all_records['basic'])

    for _, record in patient.all_records['presc'].iterrows():
        if check_presc(record, drug_name):
            presc_date = record['약품처방일']
            presc_interval = presc_date - last_adm_date

            if presc_interval > washout:
                index_dates.append(presc_date)
                
            last_adm_date = presc_date + record['투약일수'] - timedelta(1)

    return index_dates


def check_presc_conti(patient, index_dates, drug_name, max_presc_gap):
    '''
    Patient 개체를 받아서 index_dates 별로 약물 연속투여가 중단되는 날짜들(dc_dates)를 반환

    Args:
        patient (Patient): 약물 연속투여 여부를 확인할 Patient 개체
        index_dates (list of datetime): incident use가 발생한 index_date의 list
        drug_name (str): 연속투여 여부를 확인할 약물명
        max_presc_gap (int): 연속투여가 끊어지는 최대 약물 처방 간 간격 (단위: days)

    Output:
        dc_dates (list of datetime): index_date 별 약물 연속투여가 끊어진 날짜
    '''
    dc_dates = []
    max_presc_gap = timedelta(max_presc_gap)

    for idx, index_date in enumerate(index_dates):
        follow_presc_conti = False

        for _, record in patient.all_records['presc'].iterrows():
            if check_presc(record, drug_name): 
                if record['약품처방일']==index_date:
                    follow_presc_conti = True
                    last_adm_date = record['약품처방일']

                if follow_presc_conti:
                    presc_date = record['약품처방일']
                    presc_interval = presc_date - last_adm_date

                    if presc_interval >= max_presc_gap:
                        dc_dates.append(last_adm_date)
                        break

                    last_adm_date = presc_date + record['투약일수'] - timedelta(1)
        
        if len(dc_dates)<(idx + 1):
            dc_dates.append(last_adm_date)
    
    return dc_dates


def check_vict_adm(patient, victim_drug, adm_period):
    '''
    해당 patient의 index dates 이전에 victim drug을 충분히 오래 복용했는지 평가

    Args:
        patient (Patient): victim drug 복용을 평가할 대상 환자
        victim_drug (str): victim drug의 약물명
        adm_period (int): 확인하고자 하는 index date 이전의 복용 기간 길이 (단위: days)

    Return:
        is_under_victs (list of bool): 환자의 index_dates 별로 victim drug 복용 여부를 평가
    '''
    is_under_victs = []

    assert type(patient.index_dates), f'Get index dates before checking a patient under using a victim drug at the index date ({patient})'

    for index_date in patient.index_dates:    
        is_under_vict = False
        adm_period_start = index_date - timedelta(adm_period)

        for _, presc_record in patient.all_records['presc'].iterrows():
            if check_presc(presc_record, victim_drug):
                prescrp_date_start = presc_record['약품처방일']
                prescrp_date_end = prescrp_date_start + presc_record['투약일수'] - timedelta(1)

                if type(is_under_vict)==bool:
                    if prescrp_date_start <= adm_period_start:
                        if index_date <= prescrp_date_end:
                            is_under_vict = True
                            break

                        else:
                            if adm_period_start <= prescrp_date_end:
                                is_under_vict = 'need to find next records'
                                last_vict_adm = prescrp_date_end

                if type(is_under_vict)==str:                
                    if prescrp_date_start <= last_vict_adm:
                        if index_date <= prescrp_date_end:
                            is_under_vict = True
                            break
                        else:
                            last_vict_adm = prescrp_date_end

        if type(is_under_vict)==str:
            is_under_vict = False

        is_under_victs.append(is_under_vict)  

    return is_under_victs



class IncidentUserDesign():
    '''
    incident user design에 따라 연구 결과를 산출하는 객체

    Args:
        dir_supreme (str): SUPREME CDW 자료(csv파일)을 저장하는 디렉토리 주소
        study_design_filename (str): incident user design을 정리한 .xlsx 파일명
    
    Attributes:
        dir_supreme (str): SUPREME CDW 자료(csv파일)을 저장하는 디렉토리 주소
        cohort_names (list of str): design 정보를 확보한 cohort names
        study_cohort (str): 연구 결과를 산출할 대상 cohort name (cohort_names 중 선택)
        supreme_cdw (SupremeCdw): study cohort 대상의 SURPEME CDW 자료

        <아래는 개별 study cohort에서의 study design factors>
        vict_drug (str): DDI의 victim drug 약물명
        perp_drug (str): DDI의 perpetrator drug 약물명 
        cov_window (int): covariate window 길이 (days)
        washout_vict (int): victim incident use 정의에서의 washout 기간 (days)
        washout_perp(int): perpetrator incident use 정의에서의 washout 기간 (days)
        max_presc_gap_vict (int): victim 연속투여에서 최대 presciprtion gap 기간 (days)
        max_presc_gap_perp (int): perpetrator 연속투여에서 최대 presciprtion gap 기간 (days)
        otcm_follow_up (int): outcome 측정 follow-up 기간 (days)
        min_days_before_enrol (int): index date 이전의 최소 자료 확보 기간 (days)

        <아래는 PSM 및 outcome measure를 위한 진단, 처방 정보 attributes>
        cov_diag (DataFrame): 지정한 cohort에서의 covariates 진단 정보
        cov_presc (DataFrame): 지정한 cohort에서의 covariates 처방 정보
        drugs_by_group (dict): cov_presc 내 처방 정보의 세부 약물 정보
        otcm_diag (DataFrame): 지정한 cohort에서의 outcome에 진단 정보
        target_icds (dict): otcm_diag에서의 icd10_name 별 icd10_code 정보


    Methods:
        set_study_cohort: study cohort를 지정하고 이에 맞추어 design factors, supreme 자료 저장
    '''
    def __init__(self, dir_supreme, study_design_filename):
        self.dir_supreme = dir_supreme

        self.study_designs = {sheet_name: pd.read_excel(study_design_filename, 
                                                        sheet_name=sheet_name, 
                                                        engine='openpyxl')
                              for sheet_name in  ['basic', 'cov', 'otcm']}

        self.cohort_names = list(self.study_designs['basic']['cohort_name'])
        self.study_cohort = False


    def __str__(self):
        if self.study_cohort:
            return f'IncidentUserDesign(study_cohort: {self.study_cohort})'
        else:
            return f'IncidentUserDesign(study_cohort: not determined)'


    def __repr__(self):
        return self.__str__()


    def set_study_cohort(self, cohort_name):
        '''
        연구 결과를 산출할 study cohort를 지정하고, 이에 맞추어 design factors, supreme 자료 저장

        Args:
            cohort_name (str): 연구 결과를 산출할 대상 study cohort 이름
        '''

        assert cohort_name in self.cohort_names, \
            f'Please input a cohort name among {self.cohort_names}.'

        self.study_cohort = cohort_name
        self._set_design_factors(self.study_cohort)

        self.supreme_cdw = SupremeCdw(self.dir_supreme, cohort_name, drop='lab')
    

    def _set_design_factors(self, cohort_name):
        '''
        cohort name을 입력 받아서 해당 cohort에 해당하는 시험 설계 요소를 저장
        (design_setting.xlsx 파일 양식은 'trial_setting_0210.xlsx' 참고)
        
        Args: 
            cohort_name (str): 시험 설계 요소를 파악할 cohort 이름
        '''
        
        assert cohort_name in self.cohort_names, \
            f'Please input a cohort name among {self.cohort_names}.'

        #basic settings
        designs_basic = self.study_designs['basic']
        designs_basic_cohort = designs_basic[designs_basic['cohort_name']==cohort_name].reset_index()

        self.vict_drug = designs_basic_cohort['vict_name'][0]
        self.perp_drug = designs_basic_cohort['perp_name'][0]
        self.cov_window = int(designs_basic_cohort['cov_window'][0])
        self.washout_vict = int(designs_basic_cohort['washout_vict'][0])
        self.washout_perp = int(designs_basic_cohort['washout_perp'][0])
        self.max_presc_gap_vict = int(designs_basic_cohort['maximum_presc_gap_vict'][0])
        self.max_presc_gap_perp = int(designs_basic_cohort['maximum_presc_gap_perp'][0])
        self.otcm_follow_up = int(designs_basic_cohort['outcome_follow_up'][0])
        self.min_days_before_enrol = int(designs_basic_cohort['min_days_before_enrollment'][0])

        #settings for covariates
        designs_cov = self.study_designs['cov']
        designs_cov_cohort = designs_cov[[(cohort_name in c) for c in designs_cov['cohort_name']]]

        self.cov_diag = designs_cov_cohort[designs_cov_cohort['cov_type']=='diagnosis']
        self.cov_presc = designs_cov_cohort[designs_cov_cohort['cov_type']=='prescription']
        self.drugs_by_group = {drug_group: drug_names for drug_group, drug_names 
                          in zip(self.cov_presc['cov_name'], self.cov_presc['cov_for_query'])}

        #settings for outcome measures
        designs_otcm = self.study_designs['otcm']
        designs_otcm_cohort = designs_otcm[[(cohort_name in c) for c in designs_otcm['cohort_name']]]

        self.otcm_diag = designs_otcm_cohort[designs_otcm_cohort['cov_type']=='diagnosis']

        self.target_icds = {icd10_name: icd10_code for icd10_name, icd10_code 
                            in zip(self.otcm_diag['cov_name'], self.otcm_diag['cov_ICD10'])}


    def select_patients_ddi(self):
        '''
        study_cohort를 설정한 iud에서 DDI cohort 환자를 선별하여 반환

        Attributes:
            self.patients_total (list of patients): iud 내 전체 환자 정보
            self.patients_ddi (list of patients): iud 내 DDI cohort에 해당하는 환자 정보
            self.patients_not_ddi (list of patients): iud 내 DDI cohort에 속하지 않는 환자 정보
        '''
        self.patients_total = []
        self.patients_ddi = []; self.patients_not_ddi = []

        #Select patients for DDI cohort    
        print(f'Start selecting patients for DDI cohort from {self}')
        for patient in tqdm(self.supreme_cdw.set_patients()):
            #Firstly, find incidents uses for victim/perpetrator drugs and save index dates
            patient.index_dates_vict = get_index_dates_IUD(patient, 
                                                          drug_name=self.vict_drug,
                                                          washout=self.washout_vict)
            patient.index_dates_perp = get_index_dates_IUD(patient, 
                                                          drug_name=self.perp_drug,
                                                          washout=self.washout_perp)
 
            #Secondly, check continuation of victim drug prescription
            patient.dc_dates_vict = check_presc_conti(patient,
                                                      index_dates=patient.index_dates_vict,
                                                      drug_name=self.vict_drug,
                                                      max_presc_gap=self.max_presc_gap_vict)
            patient.dc_dates_perp = check_presc_conti(patient,
                                                      index_dates=patient.index_dates_perp,
                                                      drug_name=self.perp_drug,
                                                      max_presc_gap=self.max_presc_gap_perp)

            #Thirdly, identify patients for DDI cohort
            patient.index_dates_perp_ddi = find_perp_index_dates_under_vict_use(patient)
            patient.is_in_ddi, patient.index_dates_ddi, patient.dc_dates_perp_ddi = \
                find_ddi_cohort_patient(patient)

            #Lastly, classify patients
            self.patients_total.append(patient)

            if patient.is_in_ddi:
                for index_date_ddi, dc_date_perp_ddi in zip(patient.index_dates_ddi, patient.dc_dates_perp_ddi):
                    patient_pseudo = patient
                    patient_pseudo.index_date_final = index_date_ddi
                    patient_pseudo.dc_date_perp_final = dc_date_perp_ddi
                    patient_pseudo.vict_presc_days = calculate_vict_presc_days(patient_pseudo)

                    min_days = self.min_days_before_enrol
                    check_min_days = check_min_days_before_enrol(patient_pseudo, min_days)

                    if check_min_days:
                        self.patients_ddi.append(patient_pseudo)

            else:
                self.patients_not_ddi.append(patient)

        
    def get_vict_presc_days_of_ddi_cohort(self):
        '''
        DDI cohort에서 victim drug duration (vict_presc_days)를 추출

        Attribute:
            vict_presc_days_ddi (list of timedelta): DDI cohort 의 vict_presc_days 정보
        '''

        #Summarize prescription statistics of patients in DDI cohort
        self.vict_presc_days_ddi = [patient.vict_presc_days for patient in self.patients_ddi]


    def select_patients_control(self):
        '''
        study_cohort를 설정한 iud에서 control cohort 환자를 선별하여 반환

        Attributes:
            iud.patients_control (list of patients): iud 내 control cohort에 해당하는 환자 정보
        '''
        #Select patients for control cohort
        self.patients_control = []

        print('Start selecting patients for control  typing all_records in SupremeCdw.')
        for patient in tqdm(self.patients_not_ddi):
            #Firstly, determine index date of patients in the control cohort
            for index_date_vict in patient.index_dates_vict:
                for _ in range(PSUEDO_PER_PATIENT):
                    patient_pseudo = patient
                    vict_presc_days_pseudo = random.choice(self.vict_presc_days_ddi)
                    patient_pseudo.vict_presc_days = vict_presc_days_pseudo
                    patient_pseudo.index_date_final = index_date_vict + vict_presc_days_pseudo

                    min_days = self.min_days_before_enrol
                    check_min_days = check_min_days_before_enrol(patient_pseudo, min_days)

                    if check_min_days:
                        self.patients_control.append(patient_pseudo)


    def set_records_for_measurement(self):
        '''
        patients_ddi 와 patient_control에서 covariate, outcome 측정을 위한 records를 세팅

        Attributes:
            patient.all_records_cov (dictionary of Dataframe): covaraite 측정을 위한 환자별 all_records
            patient.all_records_otcm_pp (dictionary of Dataframe): outcome 측정(per-protocol)을 위한 환자별 all_records
            patient.all_records_otcm_itt (dictionary of Dataframe): outcome 측정(intetion-to-treat)을 위한 환자별 all_records
        '''               
        #set records for measuring covariates and outcomes (per-protocol and intention-to-treat)
        print('Start setting all_records for covariate and outcome measurements')
        for patient in tqdm(self.patients_ddi + self.patients_control):
            patient.all_records_cov, patient.all_records_otcm_pp, patient.all_records_otcm_itt = get_cov_otcm_records(patient, self.cov_window, self.otcm_follow_up)


def find_perp_index_dates_under_vict_use(patient):
    '''
    Patient 개체를 입력받아 victim drug prescription 중 perpetrator incident use가 발생했는지 파악하여 해당 index_dates_perp를 반환

    Arg:
        patient (Patient): index_dates_vict, index_dates_perp, dc_dates_vict을 확보한 Patient 개체

    Output:
        index_dates_perp_rev (list of datetime): victim drug prescription 중 perpetrator incident use가 발생한 index_date 리스트
    '''
    index_dates_perp_rev = []

    for idx, index_date_perp in enumerate(patient.index_dates_perp):
        for index_date_vict, dc_date in zip(patient.index_dates_vict, patient.dc_dates_vict):
            if index_date_vict<=index_date_perp<=dc_date:
                index_dates_perp_rev.append(index_date_perp)

        if len(index_dates_perp_rev)<(idx+1):
            index_dates_perp_rev.append(False)

    return index_dates_perp_rev


def find_ddi_cohort_patient(patient):
    '''
    index_dates_perp_ddi를 확보한 Patient 개체를 받아서 ddi cohort에 해당하는지 확인하고 관련 index date 정보 반환
    (현재는 첫번째 index_date_perp_ddi만 사용)

    Args:
        patient (Patient): ddi cohort 해당 여부를 판단할 Patient 개체

    Outpus:
        is_in_ddi (bool): ddi cohort에 해당하는지 여부
        index_dates_ddi (list of datetime): 확인한 ddi index date 중 첫번째 date (길이 1인 리스트로 저장)
        dc_dates_perp_ddi (list of datetime): 확인한 ddi index date에 대응하는 dc_date_perp (길이 1인 리스트로 저장)
    '''
    is_in_ddi = False
    index_dates_ddi = []
    dc_dates_perp_ddi = []

    for index_date_perp_ddi, dsct_date_perp in zip(patient.index_dates_perp_ddi, patient.dc_dates_perp):
        if index_date_perp_ddi:
            is_in_ddi = True
            index_dates_ddi.append(index_date_perp_ddi)
            dc_dates_perp_ddi.append(dsct_date_perp)
            #using only the first perpetrator index dates for identification of DDI cohort patient
            #if using all the perpetrator index dates, remove break
            break

    return is_in_ddi, index_dates_ddi, dc_dates_perp_ddi


def calculate_vict_presc_days(patient):
    '''
    final index date를  확보한 ddi cohort 내 patient를 받아 final index date에서 vict drug 복용일 계산
    
    Arg: 
        patient (Patient): index_date_final을 확보한 patient 개체

    Output:
        vict_presc_days (timedelta): victim incident use부터 final index date까지 victim drug 연속 복용 일수
    '''
    assert patient.index_date_final, f'cannot find index_date_final in a {patient}'

    vict_presc_days = timedelta(-1)

    for index_date_vict, dc_date_vict in zip(patient.index_dates_vict, patient.dc_dates_vict):
        if index_date_vict<=patient.index_date_final<=dc_date_vict:
            vict_presc_days = patient.index_date_final - index_date_vict

    assert vict_presc_days.days>=0, f'Victim drug was not continously administered at the index_date_final {vict_presc_days}'

    return vict_presc_days


def check_min_days_before_enrol(patient, min_days):
    '''
    final index date를 확정한 patient 개체에서 min_days_before_enrol 확보 여부를 확인하여 반환
    
    Args: 
        patient (Patient): final index_date가 확정된 patient 개체
        min_days (int): incident use design에서 정한 min_days_before_enrol (days)

    Output:
        check_min_days (bool): 해당 patient 개체에서의 min_days_before_enrol 확보 여부
    '''
    assert patient.index_date_final, 'Patient must obtain an attributes of index_date_final'

    min_days = timedelta(min_days)
    first_record_date = patient.all_records['basic']['최초수진일'][0]
    days_before_enrol = patient.index_date_final - first_record_date

    check_min_days = True if days_before_enrol>=min_days else False

    return check_min_days


def get_cov_otcm_records(patient, cov_period, follow_up_period):
    '''
    index_date_final을 확보한 patient 개체에서 covaraite, outcome measurment를 위한 all_records 선별

    Agr:
        patient (Patient): index_date_final 및 주요 index_date를 확보한 patient 개체
        cov_period (int): covariate masurement를 idex date 이전 covariate window 길이 (days)
        follow_up_period (int): outcome masurement를 idex date 이후 follow-up period 길이 (days)
    
    Outputs:
        all_records_cov (dictionary of DataFrame): covariate measurment를 위한 all_records
        all_records_otcm_pp (dictionary of DataFrame): outcome measurment를 위한 all_records (per-protocol, 정해진 follow-up 기간 전에 약물 연속 복용이 중단되는 경우 follow-up도 중단)
        all_records_otcm_itt (dictionary of DataFrame): outcome measurment를 위한 all_records (intention-to-treat, 약물 연속 복용 여부와 상관없이 정해진 follow-up 기간 모두 관찰)
    '''
    assert patient.index_date_final, f'cannot find index_date_final in a {patient}'
    
    try:
        patient.dc_date_perp_final
        ddi_or_control = 'ddi'
    except AttributeError:
        ddi_or_control = 'control'

    index_date_final = patient.index_date_final
    cov_period = timedelta(cov_period)
    follow_up_period = timedelta(follow_up_period)

    dc_date_vict_final = index_date_final + follow_up_period

    for index_date_vict, dc_date_vict in zip(patient.index_dates_vict, patient.dc_dates_vict):
        if index_date_vict<=index_date_final<=dc_date_vict:
            dc_date_vict_final = dc_date_vict

    if ddi_or_control=='ddi':
        dc_date_last = max(patient.dc_date_perp_final, dc_date_vict_final)
    else:
        dc_date_last = dc_date_vict_final

    censored_days = dc_date_last - index_date_final
    follow_up_period_pp = min(censored_days, follow_up_period)

    all_records_cov = patient.get_all_records_cov(cov_period.days, index_date_final)
    all_records_otcm_pp = patient.get_all_records_otcm(follow_up_period_pp.days, index_date_final)
    all_records_otcm_itt = patient.get_all_records_otcm(follow_up_period.days, index_date_final)

    return all_records_cov, all_records_otcm_pp, all_records_otcm_itt


def import_ICD10_code2name(filename):
    '''
    ICD10_cod2name 파일을 열어서 dictionary 형태로 저장
    '''
    df = pd.read_csv(filename)
    ICD10_code2name_dict = {}

    for icd10_code, icd10_name in zip(df['ICD10_Code'], df['ICD10_name']):
        ICD10_code2name_dict[icd10_code] = icd10_name

    return ICD10_code2name_dict


def get_most_common_ICDs(all_records, top_code_number, norm_level='category'):
    '''
    study design 내 전체 all_records에서 가장 빈번하게 등장한 ICD10 코드 선별하여 반환

    Args:
        all_records (dictionary of Dataframe): study design 내 전체 all_records 
         * iud.supreme_cdw.all_records 혹은 sccs.supreme_cdw.all_records
        top_code_number (int): 출현 빈도를 기준으로 선별할 상위 코드의 개수
        norm_level (int): ICD10 코드 빈도 계산 전에 수행할 정규화 수준 (2 혹은 3)

    Output:
        ICD10_codes_most_common (list of tuples): 출현 빈도를 기준으로 선별한 상위 코드 리스트
         ex. [('N18', 246200),  ('I20', 180013)]
    '''
    assert (norm_level in ['category','detail']), \
        'norm_level must be "category" or "detail"'

    diag_records = all_records['diag']
    ICD10_codes = list(diag_records['ICD10코드'])

    #normalize ICD10_codes with a given normalization_level
    if norm_level=='category':
        ICD10_codes_normed = [ICD10_code[0:3] for ICD10_code in ICD10_codes]
    else:
        ICD10_codes_normed = ICD10_codes


    #calculate the frequencies of ICD10 codes
    ICD10_codes_counter = Counter(ICD10_codes_normed)
    ICD10_codes_most_common = ICD10_codes_counter.most_common(top_code_number)

    #convert ICD10_codes_most_common into target_ICDs
    target_ICDs_most_common = {}

    for ICD10_code, _ in ICD10_codes_most_common:
        ICD10_name = ICD10_CODE2NAME_DICT[ICD10_code]
        ICD10_code_rev = (ICD10_code + 'X') if '.' in ICD10_code else (ICD10_code + '.X')
        target_ICDs_most_common[ICD10_name] = ICD10_code_rev

    return target_ICDs_most_common


def summarize_vict_presc_days(iud):
    '''
    iud 내 patients_ddi와 patients_control의 vict_presc_days 분포를 비교하여 출력

    Arg:
        iud (IncidentUserDesign): vict_presc_days를 조사할 iud
    '''
    vict_presc_days_ddi = [patient.vict_presc_days.days for patient in iud.patients_ddi]
    vict_presc_days_mean_ddi = np.mean(vict_presc_days_ddi)
    vict_presc_days_std_ddi = np.std(vict_presc_days_ddi)

    vict_presc_days_control = [patient.vict_presc_days.days for patient in iud.patients_control]
    vict_presc_days_mean_control = np.mean(vict_presc_days_control)
    vict_presc_days_std_control= np.std(vict_presc_days_control)

    print(f'vict_presc_days_mean (ddi/control): {vict_presc_days_mean_ddi:.3f}/\
        {vict_presc_days_mean_control:.3f}')
    print(f'vict_presc_days_std (ddi/control): {vict_presc_days_std_ddi:.3f}/\
        {vict_presc_days_std_control:.3f}')


#set global variables
PSUEDO_PER_PATIENT = settings.PSUEDO_PER_PATIENT
ICD10_CODE2NAME_FIELNAME = settings.ICD10_CODE2NAME_FIELNAME
ICD10_CODE2NAME_DICT = import_ICD10_code2name(ICD10_CODE2NAME_FIELNAME)

#set global variables
# DIR_BASE = settings.DIR_BASE
# STUDY_DESIGN_FILENAME = settings.STUDY_DESIGN_FILENAME

# COHORT_NAME = settings.COHORT_NAME
# ICD10_CODE_TOP_NUMBER = settings.ICD10_CODE_TOP_NUMBER
# NORM_LEVEL = settings.NORM_LEVEL
# PSUEDO_PER_PATIENT = settings.PSUEDO_PER_PATIENT

# #Incident user design for evaluating DDIs
# iud = IncidentUserDesign(dir_supreme=DIR_BASE, study_design_filename=STUDY_DESIGN_FILENAME)
# iud.set_study_cohort(COHORT_NAME)

# target_ICDs = get_most_common_ICDs(all_records=iud.supreme_cdw.all_records,
#                                    top_code_number=ICD10_CODE_TOP_NUMBER,
#                                    norm_level=NORM_LEVEL)

# iud.select_patients_ddi()
# iud.get_vict_presc_days_of_ddi_cohort()
# iud.select_patients_control()
# iud.set_records_for_measurement()

# #Compare victim drug prescription statistics of patients in ddi and control cohorts
# summarize_vict_presc_days(iud)
