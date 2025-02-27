import pandas as pd
import os

def custom_statistics(input_vcf, output_dir):
    """
    统计CustomFiltering的结果，并添加变异类型分析
    """
    # 读取文件，找到最后一个#开头的行（列名行）
    header_line = None
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_line = line.strip()
    
    if header_line is None:
        print("错误：文件中没有找到以#开头的列名行")
        return None
    
    # 使用最后一个#行作为列名来读取数据
    if header_line:
        # 处理列名中的空格，将多个空格压缩为一个
        column_names = [col.strip() for col in header_line.replace('#', '').split()]
        print("读取到的列名：", column_names)  # 调试信息
        
        # 使用空格作为分隔符读取数据，跳过标题行
        df = pd.read_csv(input_vcf, sep=r'\s+', skiprows=lambda x: x==0, names=column_names, engine='python')
        print("数据框的列名：", df.columns.tolist())  # 调试信息
        
        # 检查'VAR'列是否存在（不区分大小写）
        df.columns = [col.upper() for col in df.columns]
        print("转换为大写后的列名：", df.columns.tolist())  # 调试信息
    else:
        print("错误：找不到列名行")
        return None
    
    # 添加Types列
    def determine_type(row):
        ref = row['REF']
        var = row['VAR']
        
        # 处理包含逗号的VAR值
        if ',' in var:
            # 检查所有的VAR变异长度是否都为1
            var_variants = var.split(',')
            if len(ref) == 1 and all(len(variant) == 1 for variant in var_variants):
                return 'SBSs'
        # 处理包含星号的情况
        elif '*' in var:
            return 'SBSs'
        
        # 检查是否为SBS（长度为1的变异）
        if len(ref) == 1 and len(var) == 1:
            return 'SBSs'
        # 检查是否为插入（包含=+或+）
        elif '=+' in var or '+' in var:
            return 'INS'
        # 检查是否为缺失（包含=-或-）
        elif '=-' in var or '-' in var:
            return 'DEL'
        else:
            return 'Unknown'  # 对于无法判断的情况返回Unknown
    
    def determine_type_detail(row):
        if row['Types'] == 'SBSs':
            # 转换（Transitions）
            if row['REF'] == 'A' and row['VAR'] == 'G':
                return 'A:T>G:C', "Ts(A>G)", "Ts"
            elif row['REF'] == 'G' and row['VAR'] == 'A':
                return 'G:C>A:T', "Ts(G>A)", "Ts"
            elif row['REF'] == 'C' and row['VAR'] == 'T':
                return 'C:G>T:A', "Ts(C>T)", "Ts"
            elif row['REF'] == 'T' and row['VAR'] == 'C':
                return 'T:A>C:G', "Ts(T>C)", "Ts"
            
            # 颠换（Transversions）
            elif row['REF'] == 'A' and row['VAR'] == 'T':
                return 'A:T>T:A', "Tv(A>T)", "Tv"
            elif row['REF'] == 'A' and row['VAR'] == 'C':
                return 'A:T>C:G', "Tv(A>C)", "Tv"
            elif row['REF'] == 'A' and row['VAR'] == 'G':
                return 'A:T>G:C', "Tv(A>G)", "Tv"
            elif row['REF'] == 'T' and row['VAR'] == 'A':
                return 'T:A>A:T', "Tv(T>A)", "Tv"
            elif row['REF'] == 'T' and row['VAR'] == 'G':
                return 'T:A>G:C', "Tv(T>G)", "Tv"
            elif row['REF'] == 'T' and row['VAR'] == 'C':
                return 'T:A>C:G', "Tv(T>C)", "Tv"
            elif row['REF'] == 'C' and row['VAR'] == 'A':
                return 'C:G>A:T', "Tv(C>A)", "Tv"
            elif row['REF'] == 'C' and row['VAR'] == 'G':
                return 'C:G>G:C', "Tv(C>G)", "Tv"
            elif row['REF'] == 'G' and row['VAR'] == 'T':
                return 'G:C>T:A', "Tv(G>T)", "Tv"
            elif row['REF'] == 'G' and row['VAR'] == 'C':
                return 'G:C>C:G', "Tv(G>C)", "Tv"
            else:
                return 'SBSs', "-", "-"
            
        elif row['Types'] == 'INS':
            var = row['VAR']
            inserted_bases = ''
            if var.startswith('=+'):
                inserted_bases = var[2:]  # 获取'=+'后的序列
            elif var.startswith('+'):
                inserted_bases = var[1:]  # 获取'+'后的序列
            else:
                return 'Unknown Insertion', "-", "-"
            
            # 插入的长度就是插入序列的长度
            ins_length = len(inserted_bases)
            
            if ins_length == 1:
                return '1bp Insertion', "-", "-"
            elif 2 <= ins_length < 10:
                return '>= 2bp Insertion', "-", "-"
            else:
                return '> 10bp Insertion', "-", "-"
                
        elif row['Types'] == 'DEL':
            var = row['VAR']
            deleted_bases = ''
            if var.startswith('=-'):
                deleted_bases = var[2:]  # 获取'=-'后的序列
            elif var.startswith('-'):
                deleted_bases = var[1:]  # 获取'-'后的序列
            else:
                return 'Unknown Deletion', "-", "-"
            
            # 删除的长度就是删除序列的长度
            del_length = len(deleted_bases)
            
            if del_length == 1:
                return '1bp Deletion', "-", "-"
            elif 2 <= del_length < 10:
                return '>= 2bp Deletion', "-", "-"
            else:
                return '> 10bp Deletion', "-", "-"


    def check_ann_impact_number(row):
        try:
            ann = row.get('ANNINFO', '')
            
            if not ann:
                return 0
            
            return ann.count(',') + 1
        
        except (AttributeError, TypeError) as e:
            # 处理非字符串类型的异常情况
            print(f"Error processing ANNINFO field: {e}")
            return 0
        
    def check_Homozygous_Heterozygous(row):
        # 定义阈值常量
        HOMOZYGOUS_THRESHOLD = 75  # 因为现在是百分比值，所以使用75而不是0.75
        
        try:
            freq = row.get('FREQ')
            
            if freq is None or not isinstance(freq, (int, float)):
                raise ValueError("Invalid or missing 'FREQ' value")
            
            # 明确边界条件
            if freq >= HOMOZYGOUS_THRESHOLD:
                return "Homozygous"
            else:
                return "Heterozygous"
        
        except Exception as e:
            print(f"Error processing row: {e}")
            return "Unknown"
        

    def determine_ann_info(row):
        # 初始化一个空字典来存储解析后的信息
        parsed_info = {}
        
        # 检查row字典中是否存在'ANNINFO'键
        if 'ANNINFO' in row and row['ANNINFO']:
            parts = row['ANNINFO'].split('|')
            variant_type = parts[1]
            impact = parts[2]  
            gene_name = parts[3]
            transcript_info = parts[4] 
            variant_detail = parts[6]
            parsed_info['Variant Type'] = variant_type
            parsed_info['Impact'] = impact
            parsed_info['Gene Name'] = gene_name
            parsed_info['Transcript Info'] = transcript_info
            parsed_info['Variant Detail'] = variant_detail
            
            return parsed_info
        else:
            return {}
    

    # 在VAR列后添加Types列
    try:
        print("处理变异信息....")
        
        # 确保所有列名都是字符串类型
        df.columns = df.columns.astype(str)
        
        # 打印所有列名和它们的类型
        '''
        print("所有列名及其类型：")
        for col in df.columns:
            print(f"{col}: {type(col)}")
        '''
            
        # 尝试不同的方式找到VAR列
        var_column = None
        for col in df.columns:
            if col.strip().upper() == 'VAR':
                var_column = col
                break
                
        if var_column is None:
            raise KeyError("找不到VAR列")
            
        var_col_index = df.columns.get_loc(var_column)
        df.insert(var_col_index + 1, 'Types', df.apply(determine_type, axis=1))
        
        # 添加新的 detail 和 Ts/Tv 列
        detail_results = df.apply(determine_type_detail, axis=1)
        df.insert(var_col_index + 2, 'Detail', [x[0] if isinstance(x, tuple) else x for x in detail_results])
        df.insert(var_col_index + 3, 'Ts/Tv', [x[1] if isinstance(x, tuple) else x for x in detail_results])
        df.insert(var_col_index + 4, 'Types_detail', [x[2] if isinstance(x, tuple) else x for x in detail_results])

    except KeyError as e:
        print(f"错误：{e}")
        print("列名列表：", list(df.columns))
        return None
    
    #处理注释信息
    try:
        print("处理注释信息....")
        # 使用 'ANNINFO' 而不是 'ANN'
        ann_col_index = df.columns.get_loc("ANNINFO")
        df.insert(ann_col_index + 1, "impact_number", df.apply(check_ann_impact_number, axis=1))
    except KeyError as e:
        print(f"错误：{e}")
        print("列名列表：", list(df.columns))
        return None
    
    try:
        print("判断纯合 or 杂合？")
        # 从 COMPARISON_VALUE 列获取百分比值
        df['FREQ'] = df['COMPARISON_VALUE'].apply(lambda x: float(x.strip('%')) if isinstance(x, str) else float(x))
        df.insert(df.columns.get_loc('FREQ') + 1, "Homo_Hete", df.apply(check_Homozygous_Heterozygous, axis=1))
    except KeyError as e:
        print(f"错误：{e}")
        print("列名列表：", list(df.columns))
        return None
    except ValueError as e:
        print(f"错误：转换百分比值时出错 - {e}")
        return None
    
    # 添加新的过滤功能
    print("\n开始进行数据过滤...")
    
    # 过滤掉CHROM列以scaffold开头的行
    original_count = len(df)
    df = df[~df['CHROM'].str.startswith('scaffold', na=False)]
    scaffold_filtered = original_count - len(df)
    print(f"已过滤掉 {scaffold_filtered} 行scaffold开头的记录")
    
    # 让用户输入要过滤的样本名称
    exclude_samples = input("请输入要排除的样本名称（用英文逗号分隔，直接回车则不排除）：").strip()
    if exclude_samples:
        exclude_list = [x.strip() for x in exclude_samples.split(',')]
        sample_count_before = len(df)
        df = df[~df['SAMPLENAME'].isin(exclude_list)]
        samples_filtered = sample_count_before - len(df)
        print(f"已过滤掉 {samples_filtered} 行指定样本的记录")
    
    output_file = os.path.join(output_dir, 'Statistics_results.vcf')
    df.to_csv(output_file, sep='\t', index=False)
    
    sample_counts = df['SAMPLENAME'].value_counts()
    return sample_counts

def show_program_info():
    """显示程序信息和使用说明"""
    info = """
    ====================================================
                VCF文件变异分析处理程序
    ====================================================
    
    程序信息:
    - 版本: 1.0
    - 作者: lizhe
    - 日期: 02-26-2025
    
    主要功能:
    1. 变异类型分析（SBSs、INS、DEL）判定
    2. 转换（Ts）和颠换（Tv）统计分析
    3. 注释信息统计
    4. 纯合/杂合判断
    5. 支持过滤scaffold序列
    6. 支持按样本名称过滤
    
    使用说明:
    1. 请准备好Vascan突变检出及snpEff注释后的VCF格式的输入文件
    2. 输入文件必须包含以下列：
       CHROM, REF, VAR, COMPARISON_VALUE, ANNINFO, SAMPLENAME
    3. 程序会在当前目录生成 Statistics_results.vcf 结果文件
    4. 根据Statistics_results.vcf 结果文件判断异常样本，重新运行本程序，排除特定样
    
    ====================================================
    """
    print(info)

# 在主程序开始处添加显示信息的调用
if __name__ == "__main__":
    show_program_info()
    input_vcf = input("\n请输入待处理的文件名：")
    output_dir = "./"
    result = custom_statistics(input_vcf, output_dir)
    if result is not None:
        print("\n样本统计结果：")
        print(result)
