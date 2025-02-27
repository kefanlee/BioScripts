import pandas as pd
import os

def group_statistics(input_vcf):
    """
    按SAMPLEGROUP分组统计变异类型，并统计每个样本的详细信息
    """
    try:
        # 读取VCF文件
        df = pd.read_csv(input_vcf, sep='\t')
        
        # 确保必要的列存在
        required_columns = ['SAMPLEGROUP', 'SAMPLENAME', 'Types', 'Detail', 'Types_detail', 'Homo_Hete', "Ts/Tv"]
        if not all(col in df.columns for col in required_columns):
            raise KeyError("缺少必要的列")

        # Sheet1: 按SAMPLEGROUP分组统计
        grouped = df.groupby('SAMPLEGROUP')
        group_results = {}
        
        for group_name, group_data in grouped:
            # 计算样本数
            sample_count = len(group_data['SAMPLENAME'].unique())
            
            group_stats = {
                'Sample_Count': sample_count,
                'Total_Variants': len(group_data),
                
                # SBSs标准6种突变类型统计
                'Total_SBSs': len(group_data[group_data['Types'] == 'SBSs']),
                #转换
                'A:T>G:C': len(group_data[group_data['Detail'] == 'A:T>G:C']),
                'C:G>T:A': len(group_data[group_data['Detail'] == 'C:G>T:A']),
                #颠换
                'A:T>C:G': len(group_data[group_data['Detail'] == 'A:T>C:G']),
                'A:T>T:A': len(group_data[group_data['Detail'] == 'A:T>T:A']),
                'G:C>C:G': len(group_data[group_data['Detail'] == 'G:C>C:G']),
                'G:C>T:A': len(group_data[group_data['Detail'] == 'G:C>T:A']),
                
                # 修改后的Ts/Tv详细统计
                #转换
                'Ts(A>G)': len(group_data[group_data['Ts/Tv'] == 'Ts(A>G)']),
                'Ts(C>T)': len(group_data[group_data['Ts/Tv'] == 'Ts(C>T)']),
                #颠换
                'Tv(A>C)': len(group_data[group_data['Ts/Tv'] == 'Tv(A>C)']),
                'Tv(A>T)': len(group_data[group_data['Ts/Tv'] == 'Tv(A>T)']),
                'Tv(G>C)': len(group_data[group_data['Ts/Tv'] == 'Tv(G>C)']),
                'Tv(G>T)': len(group_data[group_data['Ts/Tv'] == 'Tv(G>T)']),
                
                'Transitions': len(group_data[group_data['Types_detail'] == 'Ts']),
                'Transversions': len(group_data[group_data['Types_detail'] == 'Tv']),
                
                # Homo_Hete统计
                'Homo': len(group_data[group_data['Homo_Hete'].str.contains('Homo', na=False)]),
                'Hete': len(group_data[group_data['Homo_Hete'].str.contains('Hete', na=False)]),

                # 保留INS详细统计
                'Total_INS': len(group_data[group_data['Types'] == 'INS']),
                '1bp_INS': len(group_data[group_data['Detail'] == '1bp Insertion']),
                '2-10bp_INS': len(group_data[group_data['Detail'] == '>= 2bp Insertion']),
                '>10bp_INS': len(group_data[group_data['Detail'] == '> 10bp Insertion']),
                
                # 保留DEL详细统计
                'Total_DEL': len(group_data[group_data['Types'] == 'DEL']),
                '1bp_DEL': len(group_data[group_data['Detail'] == '1bp Deletion']),
                '2-10bp_DEL': len(group_data[group_data['Detail'] == '>= 2bp Deletion']),
                '>10bp_DEL': len(group_data[group_data['Detail'] == '> 10bp Deletion']),
            }
            
            # 计算Ts/Tv比率
            if group_stats['Transversions'] > 0:
                group_stats['Ts/Tv_ratio'] = round(group_stats['Transitions'] / group_stats['Transversions'], 2)
            else:
                group_stats['Ts/Tv_ratio'] = 'NA'
                
            group_results[group_name] = group_stats

        # Sheet2: 按SAMPLEGROUP和SAMPLENAME分组统计
        sample_grouped = df.groupby(['SAMPLEGROUP', 'SAMPLENAME'])
        sample_results = []
        
        for (group_name, sample_name), sample_data in sample_grouped:
            sample_stats = {
                'SAMPLEGROUP': group_name,
                'SAMPLENAME': sample_name,
                'Total_Variants': len(sample_data),
                'Total_SBSs': len(sample_data[sample_data['Types'] == 'SBSs']),
                
                # 标准突变类型统计
                'A:T>G:C': len(sample_data[sample_data['Detail'] == 'A:T>G:C']),
                'C:G>T:A': len(sample_data[sample_data['Detail'] == 'C:G>T:A']),
                'A:T>C:G': len(sample_data[sample_data['Detail'] == 'A:T>C:G']),
                'A:T>T:A': len(sample_data[sample_data['Detail'] == 'A:T>T:A']),
                'G:C>C:G': len(sample_data[sample_data['Detail'] == 'G:C>C:G']),
                'G:C>T:A': len(sample_data[sample_data['Detail'] == 'G:C>T:A']),
                
                # Ts/Tv详细统计
                'Ts(A>G)': len(sample_data[sample_data['Ts/Tv'] == 'Ts(A>G)']),
                'Ts(C>T)': len(sample_data[sample_data['Ts/Tv'] == 'Ts(C>T)']),
                'Tv(A>C)': len(sample_data[sample_data['Ts/Tv'] == 'Tv(A>C)']),
                'Tv(A>T)': len(sample_data[sample_data['Ts/Tv'] == 'Tv(A>T)']),
                'Tv(G>C)': len(sample_data[sample_data['Ts/Tv'] == 'Tv(G>C)']),
                'Tv(G>T)': len(sample_data[sample_data['Ts/Tv'] == 'Tv(G>T)']),
                
                'Transitions': len(sample_data[sample_data['Types_detail'] == 'Ts']),
                'Transversions': len(sample_data[sample_data['Types_detail'] == 'Tv']),
                
                'Homo': len(sample_data[sample_data['Homo_Hete'].str.contains('Homo', na=False)]),
                'Hete': len(sample_data[sample_data['Homo_Hete'].str.contains('Hete', na=False)]),
                
                # INS/DEL详细统计
                'Total_INS': len(sample_data[sample_data['Types'] == 'INS']),
                '1bp_INS': len(sample_data[sample_data['Detail'] == '1bp Insertion']),
                '2-10bp_INS': len(sample_data[sample_data['Detail'] == '>= 2bp Insertion']),
                '>10bp_INS': len(sample_data[sample_data['Detail'] == '> 10bp Insertion']),
                'Total_DEL': len(sample_data[sample_data['Types'] == 'DEL']),
                '1bp_DEL': len(sample_data[sample_data['Detail'] == '1bp Deletion']),
                '2-10bp_DEL': len(sample_data[sample_data['Detail'] == '>= 2bp Deletion']),
                '>10bp_DEL': len(sample_data[sample_data['Detail'] == '> 10bp Deletion'])
            }
            
            # 计算单个样本的Ts/Tv比率
            if sample_stats['Transversions'] > 0:
                sample_stats['Ts/Tv_ratio'] = round(sample_stats['Transitions'] / sample_stats['Transversions'], 2)
            else:
                sample_stats['Ts/Tv_ratio'] = 'NA'
                
            sample_results.append(sample_stats)
        
        # 转换为DataFrame
        group_df = pd.DataFrame.from_dict(group_results, orient='index')
        sample_df = pd.DataFrame(sample_results)
        
        # 保存到Excel文件，创建两个sheet
        output_file = 'group_statistics_results.xlsx'
        with pd.ExcelWriter(output_file) as writer:
            group_df.to_excel(writer, sheet_name='Group_Statistics')
            sample_df.to_excel(writer, sheet_name='Sample_Statistics', index=False)
            
        print(f"结果已保存到 {output_file}")
        print("\n分组统计结果：")
        print(group_df)
        print("\n样本统计结果：")
        print(sample_df)
        
        return group_df, sample_df
        
    except Exception as e:
        print(f"处理过程中出现错误：{str(e)}")
        return None, None

def show_program_info():
    """显示程序信息和使用说明"""
    info = """
    ====================================================
                VCF文件分组统计分析程序
    ====================================================
    
    程序信息:
    - 版本: 1.0
    - 作者: lizhe
    - 日期: 02-26-2025
    
    主要功能:
    1. 按样本组(SAMPLEGROUP)进行分组统计
    2. 统计每个样本的详细变异信息
    3. SNP的6种标准突变类型统计
    4. 转换（Ts）和颠换（Tv）比率计算
    5. 纯合/杂合突变数量统计
    6. 插入(INS)和缺失(DEL)的详细分类统计
    
    使用说明:
    1. 请准备好第一步生成的 Statistics_results.vcf 文件
    2. 输入文件必须包含以下列：
       SAMPLEGROUP, SAMPLENAME, Types, Detail, Types_detail, 
       Homo_Hete, Ts/Tv
    3. 程序会生成 group_statistics_results.xlsx 结果文件，包含：
       - Group_Statistics：样本组整体统计结果
       - Sample_Statistics：单个样本详细统计结果

    特别说明：程序生成的group_statistics_results.xlsx 结果文件可用于 `R` 绘图统计分析。
    
    ====================================================
    """
    print(info)

if __name__ == "__main__":
    show_program_info()
    input_file = input("请输入待处理的文件名(默认：Statistics_results.vcf，回车即可)：")
    if input_file == "":
        input_file = "Statistics_results.vcf"
    if os.path.exists(input_file):
        group_statistics(input_file)
    else:
        print("当前目录下未找到 Statistics_results.vcf 文件") 