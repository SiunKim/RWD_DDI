import pandas as pd

from datetime import timedelta
from collections import defaultdict


DATE_COLUMNS = ['수진(진료)일', '약품처방일', '진단일자', '검사시행일']
RECORD_TYPES_ABBR = {'basic': 'basic', 'exam': 'examination', 'presc': 'prescription', 'diag': 'diagnosis', 'lab': 'labtest'}
# COLUMNS = {'exam': ['원무접수ID', '환자번호', '성별', '생년월일', '수진(진료)일', '수진(퇴원포함)진료과', '진료구분(초/재진/신환)', '입원일', '퇴원일']} 


def check_presc(record, drug_name):
    def check_presc_indi(record, drug_name):
        parsed_drug_names = record['약품명(성분명)'].lower().split('/')
        is_presced = any([parsed_drug_name.startswith(drug_name.lower()) 
                          for parsed_drug_name in parsed_drug_names])

        return is_presced

    '''
    개별 record line에서 약물 복용 여부를 확인해서 True, False 반환

    Args:
        record (DataFrame): prescription 정보를 담은 record DataFrame 한 행
        drug_name (str): 복용 여부를 확인할 약물명 (drug    group의 경우는 '<group>'으로 시작)
    
    Return:
        is_presced (bool): 복용 여부
    '''
    if drug_name.startswith('<group>'):
        group_name = drug_name.replace('<group>', '')

        assert group_name in DRUGS_BY_GROUP.keys(), \
        f'<group>group_name(group_name{group_name}) is not defined in DRUGS_BY_GROUP '

        is_presced = any([check_presc_indi(record, drug_ind)
                          for drug_ind in DRUGS_BY_GROUP[group_name].split(',')])

        return is_presced

    else:
        is_presced = check_presc_indi(record, drug_name)

        return is_presced


def records_between(records, date_start, date_end):
    '''
    개별 records를 받아서 주어진 기간에 맞추어 filtering하고 반환
    '''
    date_cols = [col for col in records.columns if col in DATE_COLUMNS]

    assert len(date_cols)==1, f'Please check the column names of records ({records.columns})'
    date_col = date_cols[0]

    is_between = pd.Series([(date_start <= date <= date_end) 
                    for date in records[date_col]], dtype=bool)

    records_between = records[is_between].reset_index(drop=True)

    return records_between


class Patient():
    '''
    환자 개인의 CDW 정보를 저장하고 때에 맞게 전달합니다.

    Args:
        ptnt_id (int): 해당 환자의 환자 번호
        all_records_ind (dictionary of DataFrame): 환자 개인의 CDW 자료

    Attributes:
        id (int): CDW 내 '환자번호'
        all_records (dictionary of DataFrame): 환자 개인의 CDW 자료
        record_period (tuple of datetiem): 
        index_dates (list of datetime): 

    Methods:
        filter_all_records: index_date 전후로 all_records 자료 filtering
        get_all_records_cov: covariate 평가를 위해 index_date 이전 all_records 자료 filtering하고 all_records_cov에 저장
        get_all_records_otcm: outcome 평가를 위해 index_date 이후 all_records 자료 filtering하고 all_records_otcm에 저장
    '''
    def __init__(self, ptnt_id, all_records_ind):
        self.id = ptnt_id
        self.all_records = all_records_ind


    def __str__(self):
        return f"Patient(id:{self.id})"


    def __repr__(self):
        return f"Patient(id:{self.id})"


    def filter_all_records(self, period_length, index_date, previous=True):
        '''
        index date 전후로 all_records 자료 filtering

        Args:
            period_length (int): index date 전후로 CDW 자료를 filtering할 기간 (단위: days)
            index_date (datetime): CDW 자료를 filtering할 기준 index date
            previous (bool): index date 이전/이후 (True/False) 날짜에 해당하는 자료를 반환

        Return:
            all_records_selected (dictionary of DataFrame): 지정 조건에 맞춰서 filtering한 새로운 all_records
        '''
        period_length = timedelta(days=period_length)

        all_records_filtered = {}

        for record_type, record in self.all_records.items():
            if previous:
                date_start = index_date - period_length
                date_end = index_date
            else:
                date_start = index_date
                date_end = index_date + period_length

            all_records_filtered[record_type] = records_between(record, date_start, date_end)

        return all_records_filtered


    def get_all_records_cov(self, cov_period, index_date):
        '''
        covariate 평가를 위해 index_date 이전 all_records 자료 filtering하고 all_records_cov에 저장
        '''
        return self.filter_all_records(cov_period, index_date)


    def get_all_records_otcm(self, follow_up_period, index_date):
        '''
        outcome 평가를 위해 index_date 이후 all_records 자료 filtering하고 all_records_outcome 저장
        '''
        return self.filter_all_records(follow_up_period, index_date, previous=False)