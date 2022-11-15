import pandas as pd
import numpy as np

from collections import defaultdict
from datetime import timedelta, datetime

from patient import check_presc, records_between
from supreme_cdw import SupremeCdw
from incident_user_design import calculate_vict_presc_days



def get_index_dates(patient, drug_name, washout):
    '''
    Patient 개체를 받아서 약물의 index_dates를 반환

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



class SelfControlledCaseSeries():
    '''
    self-controlled case series design에 따라 연구 결과를 산출하는 객체

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
        cov_window (int): index date 이전 covariate 측정을 위해 추출하는 records 기간 (days) - 시험 설계 상으로 SCCS에서는 필요하지 않으나 ICU와 환자 특성 비교를 위해 정의
        washout_vict (int): victim incident use 정의에서의 washout 기간 (days)
        max_presc_gap_vict (int): victim 연속투여에서 최대 presciprtion gap 기간 (days)
        max_presc_gap_perp (int): perpetrator 연속투여에서 최대 presciprtion gap 기간 (days)
        min_days_before_enrol (int): index date 이전의 최소 자료 확보 기간 (days)
        induction_window (int): perpetrator 투여 시작 이후 perp-exposed 기간 시작까지 간격 (days)
        post_exposure_followup (int): perpetrator 연속투여 중단 이후 perp-exposed 추척 기간 (days)
        post_exposure_washout (int): post-exposure follow-up 기간 종료 이후 washout 기간 (days)

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

        self.supreme_cdw = SupremeCdw(dir_base, cohort_name, drop='lab')
    

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
        self.max_presc_gap_vict = int(designs_basic_cohort['maximum_presc_gap_vict'][0])
        self.max_presc_gap_perp = int(designs_basic_cohort['maximum_presc_gap_perp'][0])
        self.min_days_before_enrol = int(designs_basic_cohort['min_days_before_enrollment'][0])
        self.induction_window = int(designs_basic_cohort['induction_window'][0])
        self.post_exposure_followup = int(designs_basic_cohort['post_exposure_followup'][0])
        self.post_exposure_washout = int(designs_basic_cohort['post_exposure_washout'][0])

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


def check_min_days_before_enrol(patient, min_days):
    '''
    final index date를 확정한 patient 개체에서 min_days_before_enrol 확보 여부를 확인하여 반환
    
    Args: 
        patient (Patient): final index_date가 확정된 patient 개체
        min_days (int): incident use design에서 정한 min_days_before_enrol (days)

    Output:
        check_min_days (bool): 해당 patient 개체에서의 min_days_before_enrol 확보 여부
    '''
    assert(patient.index_date_final, 'Patient must obtain an attributes of index_date_final')

    min_days = timedelta(min_days)
    first_record_date = patient.all_records['basic']['최초수진일'][0]
    days_before_enrol = patient.index_date_final - first_record_date

    check_min_days = True if days_before_enrol>=min_days else False

    return check_min_days


def check_and_correct_followups_starts_ends(followups_starts_ends):
    '''
    followups_starts_ends에서 start, end date의 수가 동일하고 start date가 end date보다 빠른지 확인하고 올바른 형식이 아닌 경우 수정하여 반환  

    Outputs:
        followups_starts_ends_rev (dictionary, list of datetime): 올바른 형식이 아닌 경우 수정한 followups_starts_ends
        checked (bool): unexposed, exposed에서 모두 start, end date의 수가 동일한지 여부
        corrected (bool): unexposed, exposed에서 start date가 end date보다 빠르지 않은 경우가 존재하는지 여부 (해당 start, end dates pair를 제거하여 followups_starts_ends_rev 반환)
    '''
    followups_starts_ends_rev = {'perp_unexposed_start': [], 'perp_unexposed_end': [],
                             'perp_exposed_start': [], 'perp_exposed_end': []}
    checked = True
    corrected = False

    unexposed_starts = followups_starts_ends['perp_unexposed_start']
    unexposed_ends = followups_starts_ends['perp_unexposed_end']
    exposed_starts = followups_starts_ends['perp_exposed_start']
    exposed_ends = followups_starts_ends['perp_exposed_end']

    if len(unexposed_starts)!=len(unexposed_ends):
        checked = False
    else:
        for unexposed_start, unexposed_end in zip(unexposed_starts, unexposed_ends):
            if unexposed_start<unexposed_end:
                followups_starts_ends_rev['perp_unexposed_start'].append(unexposed_start)
                followups_starts_ends_rev['perp_unexposed_end'].append(unexposed_end)
            else:
                corrected = True
                
    if len(exposed_starts)!=len(exposed_ends):
        checked = False
    else:
        for exposed_start, exposed_end in zip(exposed_starts, exposed_ends):
            if exposed_start<exposed_end:
                followups_starts_ends_rev['perp_exposed_start'].append(exposed_start)
                followups_starts_ends_rev['perp_exposed_end'].append(exposed_end)
            else:
                corrected = True

    return checked, corrected, followups_starts_ends_rev


def arrange_presc_periods(presc_periods, max_presc_gap=False):
    '''
    presc_periods를 받아 중복되는 기간을 정리하여 다시 반환

    Args:
        presc_periods (dictionary of list of datatime): 정리할 presc_periods dictionary
         - keys: 'period_start', 'period_end'
        max_presc_gap_perp (int, optional): 연속투여가 끊어지는 최대 약물 처방 간 간격 (단위: days)

    Outpus:
        presc_periods_rev (dictionary of list of datatime): 중복 기간을 정리한 presc_periods
    '''
    period_starts = presc_periods['period_start']
    period_ends = presc_periods['period_end']

    if period_starts:
        last_period_start = period_starts[0]
        last_period_end = period_ends[0]

        max_presc_gap = timedelta(max_presc_gap) if max_presc_gap else False

        presc_periods_rev = defaultdict(list)
        presc_periods_rev['period_start'].append(last_period_start)
        presc_periods_rev['period_end'].append(last_period_end)
        
        for period_start, period_end in zip(period_starts[1:], period_ends[1:]):
            if max_presc_gap:
                is_overlapped = bool(last_period_end + max_presc_gap + timedelta(1) >= period_start)
            else:
                is_overlapped = bool(last_period_end + timedelta(1) >= period_start)

            if is_overlapped:
                presc_periods_rev['period_end'].pop()
                presc_periods_rev['period_end'].append(period_end)

            else:
                presc_periods_rev['period_start'].append(period_start)
                presc_periods_rev['period_end'].append(period_end)

            last_period_start = period_start
            last_period_end = period_end

    else:
        presc_periods_rev = presc_periods

    return presc_periods_rev


def extract_presc_periods(presc_records, drug_name, max_presc_gap=False):
    '''
    presc_records에서 drug_name에 해당하는 약물의 복용일을 추출하여 반환
    Args:
        presc_records (Dataframe): 약물 복용일을 확인할 records Dataframe
        drug_name (str): 복용일을 확인할 약물명
        max_presc_gap_perp (int, optional): 연속투여가 끊어지는 최대 약물 처방 간 간격 (단위: days)

    Output:
        presc_periods (dictionary of list of datetime): 약물 복용일의 시작일과 종료일을 각각 list로 저장하여 dictionary 구성
        keys: 'period_start', 'period_end'
    '''
    presc_periods = defaultdict(list)

    for _, record in presc_records.iterrows():
        if check_presc(record, drug_name):
            presc_date = record['약품처방일']
            presc_periods['period_start'].append(presc_date)
            presc_periods['period_end'].append(presc_date + record['투약일수'] - timedelta(1))

    # presc_periods_rev = arrange_presc_periods(presc_periods, max_presc_gap)

    return presc_periods


def transform_presc_periods_to_periodindex(presc_periods):
    '''
    presc_periods를 pandas PeriodIndex로 변환
    '''
    period_starts = presc_periods['period_start']
    period_ends = presc_periods['period_end']

    if period_starts:
        periodindex_total = pd.period_range(period_starts[0], period_ends[0])

        if len(period_starts)>1:
            for period_start, period_end in zip(period_starts[1:], period_ends[1:]):
                periodindex = pd.period_range(period_start, period_end)
                periodindex_total = periodindex_total.union(periodindex)

    else:
        periodindex_void = pd.period_range(datetime.today(), datetime.today() - timedelta(1))
        periodindex_total = periodindex_void

    return periodindex_total

    
def transform_periodindex_to_presc_periods(periodindex):
    # def period_to_datetime(period):
    #     return datetime.fromtimestamp(period.to_timestamp().timestamp())
    '''
    pandas PeriodIndex를 presc_periods로 변환
    '''
    presc_periods = defaultdict(list)

    if bool(list(periodindex)):
        last_datetime = periodindex[0].to_timestamp()
        presc_periods['period_start'].append(last_datetime)

        for period in periodindex[1:-1]:
            current_datetime = period.to_timestamp()

            #when consequent periodindexes are discontinued
            if last_datetime + timedelta(1)!=current_datetime:
                presc_periods['period_end'].append(last_datetime)
                presc_periods['period_start'].append(current_datetime)

            last_datetime = current_datetime

        end_datetime = periodindex[-1].to_timestamp()
        presc_periods['period_end'].append(end_datetime)

    else:        
        presc_periods['period_start'] = []
        presc_periods['period_end'] = []

    return presc_periods


def get_union_of_presc_periods(presc_periods1, presc_periods2):
    '''
    두 presc_periods을 합쳐서 하나의 presc_period로 반환
    '''
    periodindex1 = transform_presc_periods_to_periodindex(presc_periods1)
    periodindex2 = transform_presc_periods_to_periodindex(presc_periods2)

    periodindex_union = periodindex1.union(periodindex2)
    presc_periods_union = transform_periodindex_to_presc_periods(periodindex_union)

    return presc_periods_union

    
def get_difference_of_presc_periods(presc_periods1, presc_periods2):
    '''
    두 presc_periods를 빼서 하나의 presc_period로 반환 (presc_periods1 - presc_periods2)
    '''
    periodindex1 = transform_presc_periods_to_periodindex(presc_periods1)
    periodindex2 = transform_presc_periods_to_periodindex(presc_periods2)

    periodindex_difference = periodindex1.difference(periodindex2)
    presc_periods_difference = transform_periodindex_to_presc_periods(periodindex_difference)

    return presc_periods_difference


def get_perp_exposed_from_perp_presc_periods(perp_presc_periods, 
                                             induction_window,
                                             post_exposure_followup):
    '''
    perp_presc_periods로부터 perp_exposed_periods를 계산하여 반환

    Args:
        perp_presc_periods (dictionary of list of datetime): perp drug 복용 기간 (presc_periods)
        induction_window (int): perpetrator 투여 시작 이후 perp-exposed 기간 시작까지 간격 (days)
        post_exposure_followup (int): perpetrator 연속투여 중단 이후 perp-exposed 추척 기간 (days)

    Output:
        perp_exposed_periods (dictionary of list of datetime): DDI 평가 시 perp-exposued period로 분류되는 기간 (presc_periods)
    '''    
    period_starts = perp_presc_periods['period_start']
    period_ends = perp_presc_periods['period_end']
    induction_window = timedelta(induction_window)
    post_exposure_followup = timedelta(post_exposure_followup)

    perp_exposed_periods = defaultdict(list)

    for period_start, period_end in zip(period_starts, period_ends):
        perp_exposed_periods['period_start'].append(period_start + induction_window)
        perp_exposed_periods['period_end'].append(period_end + post_exposure_followup)

    perp_exposed_periods_rev = arrange_presc_periods(perp_exposed_periods)

    return perp_exposed_periods_rev

    
def get_washout_from_perp_presc_periods(perp_presc_periods,
                                        post_exposure_followup,
                                        post_exposure_washout):
    '''
    perp_presc_periods로부터 wahout_periods를 계산하여 반환

    Args:
        perp_presc_periods (dictionary of list of datetime): perp drug 복용 기간 (presc_periods)
        post_exposure_followup (int): perpetrator 연속투여 중단 이후 perp-exposed 추척 기간 (days)
        post_exposure_washout (int): post-exposure follow-up 기간 종료 이후 washout 기간 (days)

    Output:
        wahout_periods를 (dictionary of list of datetime): DDI 평가 시 post-exposure washout period로 분류되는 기간 (presc_periods)
    '''    
    period_starts = perp_presc_periods['period_start']
    period_ends = perp_presc_periods['period_end']
    post_exposure_followup = timedelta(post_exposure_followup)
    post_exposure_washout = timedelta(post_exposure_washout)

    wahout_periods = defaultdict(list)

    for _, period_end in zip(period_starts, period_ends):
        washout_period_start = period_end + post_exposure_followup + timedelta(1)
        washout_period_end = washout_period_start + post_exposure_washout

        wahout_periods['period_start'].append(washout_period_start)
        wahout_periods['period_end'].append(washout_period_end)

    wahout_periods_rev = arrange_presc_periods(wahout_periods)

    return wahout_periods_rev


def sort_out_total_followup_period(patient,
                                   perp_drug, 
                                   induction_window,
                                   max_presc_gap_perp,
                                   post_exposure_followup,
                                   post_exposure_washout):
    '''
    index_date_final, dc_date_vict_final를 확보한 Patient 개체를 받아서 해당 follow-up period를 perp-unexposed, perp-exposed, and post-exposure washout periods별로 분류

    Args:
        patient (Patient): follow-up period를 분류할 Patient 개체
        perp_drug (str): DDI의 perpetrator drug 약물명 
        induction_window (int): perpetrator 투여 시작 이후 perp-exposed 기간 시작까지 간격 (days)
        max_presc_gap_perp (int): 연속투여가 끊어지는 최대 약물 처방 간 간격 (단위: days)
        post_exposure_followup (int): perpetrator 연속투여 중단 이후 perp-exposed 추척 기간 (days)
        post_exposure_washout (int): post-exposure follow-up 기간 종료 이후 washout 기간 (days)

    Output:
        followups_starts_ends (dictionary, list of datetime): perp-unexposed, perp-exposed, post-exposure washout period 별로 시작일/종료일(datetime) list로 저장
        keys: 'perp_unexposed_start', 'perp_unexposed_end', 'perp_exposed_start', 'perp_exposed_end'
    '''
    followup_start = patient.index_date_final
    followup_end = patient.dc_date_vict_final
    presc_records = records_between(patient.all_records['presc'], followup_start, followup_end)
        
    # set total_followup_preiods
    total_followup_period = {'period_start': [followup_start], 'period_end': [followup_end]}

    # extact perp drug presciption periods
    perp_presc_periods = extract_presc_periods(presc_records, perp_drug, max_presc_gap_perp)
    # get perp exposed periods from perp drug presciption periods (considering induction window, post-exposure followup)
    perp_exposed_periods = get_perp_exposed_from_perp_presc_periods(perp_presc_periods,
                                                                    induction_window,
                                                                    post_exposure_followup)
    # get perp washout periods from perp drug presciption periods (considering induction window, post-exposure followup)
    perp_wahout_periods = get_washout_from_perp_presc_periods(perp_presc_periods, 
                                                              post_exposure_followup,
                                                              post_exposure_washout)


    not_perp_exposed_periods1 = get_difference_of_presc_periods(total_followup_period,
                                                             perp_presc_periods)
    not_perp_exposed_periods2 = get_difference_of_presc_periods(not_perp_exposed_periods1,
                                                             perp_exposed_periods)
    perp_unexposed_periods = get_difference_of_presc_periods(not_perp_exposed_periods2,
                                                             perp_wahout_periods)

    return perp_unexposed_periods, perp_exposed_periods


def summarize_patient_ddi(patient_ddi):
    '''
    patient_ddi 환자의 perp_unexposed_periods, perp_exposed_periods 길이를 계산해 반환
    '''
    unexposed_periods_starts = patient_ddi.perp_unexposed_periods['period_start']
    unexposed_periods_ends = patient_ddi.perp_unexposed_periods['period_end']
    exposed_periods_starts = patient_ddi.perp_exposed_periods['period_start']
    exposed_periods_ends = patient_ddi.perp_exposed_periods['period_end']

    unexposed_periods_length = timedelta(0)
    exposed_periods_length = timedelta(0)

    for period_start, period_end in zip(unexposed_periods_starts, unexposed_periods_ends):
        unexposed_periods_length += (period_end - period_start)
    for period_start, period_end in zip(exposed_periods_starts, exposed_periods_ends):
        exposed_periods_length += (period_end - period_start)

    return unexposed_periods_length, exposed_periods_length


def test_patient_ddi(sccs, patient_ddi_index):
    '''
    개별 patient_ddi에서 attributes 확인 및 출력
    '''
    patient_ddi = sccs.patients_ddi[patient_ddi_index]

    print(f'perp: {sccs.perp_drug}, vict: {sccs.vict_drug}')
    print(f'induction_window: {sccs.induction_window}, post_exposure_followup: {sccs.post_exposure_followup}, post_exposure_washout: {sccs.post_exposure_washout}')
    
    print(f'perp_unexposed_periods: {patient_ddi.perp_unexposed_periods}')
    print(f'perp_exposed_periods: {patient_ddi.perp_exposed_periods}')
    print(f'index_dates_vict: {patient_ddi.index_dates_vict}')
    print(f'index_date_final: {patient_ddi.index_date_final}')
    print(f'dc_date_vict_final: {patient_ddi.dc_date_vict_final}')

    presc_records = records_between(patient_ddi.all_records['presc'],
                                    patient_ddi.index_date_final, 
                                    patient_ddi.dc_date_vict_final)

    display(presc_records)
    

#sample codes
dir_base = 'D:/Google_Drive/[1] CCADD/[1] Personal Research/2021_DDI_RWD\DATASET'
filename = 'trial_setting_0210.xlsx'

sccs = SelfControlledCaseSeries(dir_supreme=dir_base, study_design_filename=filename)
sccs.cohort_names

sccs.set_study_cohort('clo_ome')

sccs.patients_total = []
sccs.patients_ddi = []; sccs.patients_no_ddi = []

#Select patients for the self-controlled case series design
for idx, patient in enumerate(sccs.supreme_cdw.set_patients()):
    #Firstly, find incidents uses for victim drug and save index dates
    patient.index_dates_vict = get_index_dates(patient, 
                                               drug_name=sccs.vict_drug,
                                               washout=sccs.washout_vict)

    #Secondly, check continuation of prescription of victim drugs
    patient.dc_dates_vict = check_presc_conti(patient,
                                              index_dates=patient.index_dates_vict,
                                              drug_name=sccs.vict_drug,
                                              max_presc_gap=sccs.max_presc_gap_vict)

    #Thirdly, sort out total outcome follow-up period into perp-unexposed, perp-exposed, and post-exposure washout periods in each patient
    if patient.index_dates_vict:
        for index_date_vict, dc_date_vict in zip(patient.index_dates_vict, patient.dc_dates_vict):
            patient_pseudo = patient
            patient_pseudo.index_date_final = index_date_vict
            patient_pseudo.dc_date_vict_final = dc_date_vict
            patient_pseudo.perp_unexposed_periods, patient_pseudo.perp_exposed_periods = \
                sort_out_total_followup_period(patient_pseudo,
                                               perp_drug=sccs.perp_drug, 
                                               induction_window=sccs.induction_window,
                                               max_presc_gap_perp=sccs.max_presc_gap_perp,
                                               post_exposure_followup=sccs.post_exposure_followup,
                                               post_exposure_washout=sccs.post_exposure_washout)

        min_days = sccs.min_days_before_enrol
        check_min_days = check_min_days_before_enrol(patient_pseudo, min_days)

        find_unexposed = bool(patient_pseudo.perp_unexposed_periods['period_start'])
        find_exposed = bool(patient_pseudo.perp_exposed_periods['period_start'])

        patient_pseudo.unexposed_periods_length, patient_pseudo.exposed_periods_length = \
            summarize_patient_ddi(patient_pseudo)
        patient_pseudo.vict_presc_days = calculate_vict_presc_days(patient_pseudo)

        if check_min_days:
            if find_unexposed & find_exposed:
                sccs.patients_ddi.append(patient_pseudo)
            else:
                sccs.patients_no_ddi.append(patient_pseudo)

    #Lastly, summarize follow-up period statistics
    sccs.patients_total.append(patient)


def filter_all_records_by_presc_peroids(patient, presc_periods):
    '''
    patient의 all_records 중 presc_periods에 해당하는 날짜의 record만을 추출하여 반환
    '''
    period_starts = presc_periods['period_start']
    period_ends = presc_periods['period_end']

    all_records = patient.all_records

    all_records_rev = defaultdict(pd.DataFrame)

    for records_type, records in all_records.items():
        for period_start, period_end in zip(period_starts, period_ends):
            records_tmp = records_between(records, period_start, period_end)
            all_records_rev[records_type] = pd.concat([all_records_rev[records_type], records_tmp])

    return all_records_rev


def get_all_records_perp_unexposed(patient):
    return filter_all_records_by_presc_peroids(patient, patient.perp_unexposed_periods)


def get_all_records_perp_exposed(patient):
    return filter_all_records_by_presc_peroids(patient, patient.perp_exposed_periods)


#Filter all_records for covariates and outcomes measurement
COV_WINDOW = 180 #days
for patient_ddi in sccs.patients_ddi:
    patient_ddi.all_records_perp_unexposed = get_all_records_perp_unexposed(patient_ddi)
    patient_ddi.all_records_perp_exposed = get_all_records_perp_exposed(patient_ddi)

    patient_ddi.all_records_cov = patient_ddi.get_all_records_cov(COV_WINDOW, 
                                                                  patient_ddi.index_date_final)
