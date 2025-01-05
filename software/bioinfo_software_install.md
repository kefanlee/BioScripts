> 本文记录了**基因组组装、全基因组重测序、转录组**数据分析等生信分析中的环境安装

# 基因组分析

## Structure Variant

### Delly

```bash
github 下载编译包-创建软连接即可
```

### Lumpy

```bash
conda create -n lumpy python=2.7 
 
conda install -c bioconda lumpy-sv=0.3.1
```

> 创建软连接
```
lumpy
```

### Manta

```bash
conda create -n manta python=2.7

conda install -c bioconda manta=1.6.0
```

> 创建软连接
```
configManta.py
runMantaWorkflowDemo.py
```