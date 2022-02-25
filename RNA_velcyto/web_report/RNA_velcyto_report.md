## 10X genomics单细胞转录组RNA速度分析结题报告

### 摘要

使用Velocyto及scVelo对{N}个单细胞样本进行RNA速度分析，预测细胞的未来状态。

### 1.背景介绍

单细胞测序需要通过破坏细胞以获取单个细胞中的基因活性信息，这使得研究细胞动态过程以及了解细胞决策具有挑战性。La Mano等人引入RNA速度(RNA velocity)，即利用新转录的，未发生剪接的前体mRNA与成熟的，发生剪接的mRNA，可在单细胞转录组数据中表征细胞定向的动态信息（通过检测内含子的存在与否）。“RNA速度”不仅可以表征基因的活性，而且可以描述基因在单个细胞中的变化，作者验证了在神经嵴谱系中的准确性，并在多个已发表的数据集和技术平台上验证了其准确性，开辟了研究细胞分化的新方法。

Bergen等基于似然的动力学模型开发的scVelo进一步将“RNA-速度”推广到与细胞瞬时状态相关的各种系统（如细胞发育以及对扰动的应激）。该动力模型解决了拼接动力学每个基因的全部动力（最强大而计算最昂贵）。使用scVelo可以推断包括转录，剪接和降解在内的”基因特异性速度(gene-specific rates)，并且可还原细胞过程中的"潜在时间(latent-time)"，而潜伏时间近为细胞分化的真实时间。更重要的是，scVelo可以识别调控机制并系统性的推定驱动基因。

该分析研究意义主要在于：

1. 根据检测细胞的类型推断细胞分化的轨迹

2) 预测单个细胞的变化方向，得到细胞间的转变过程
   
3) 对以往数据进行挖掘，采用初生（未剪接的）和成熟（剪接的）mRNA的相对丰度用来估计基因剪接和降解的速率
4) 未稳定的细胞向处于稳泰平衡的细胞分化，揭示动力学稳态研究过程

![image-20211217161437097](C:\Users\jmzha\AppData\Roaming\Typora\typora-user-images\image-20211217161437097.png)

- a.转录动力学建模捕捉未剪接前 mRNA 的转录激活和抑制（“开”和“关”阶段）这些转录本转化为成熟的mRNA并最终降解；b.当转录激活和抑制的转录阶段分别持续足够长的时间时，分别达到活跃转录和非活跃沉默的稳态。然而，特别是在瞬态细胞群中，通常不会达到稳定状态，例如，转录激活可能在 mRNA 水平饱和之前终止，表现出"提早转换”行为。c. 基于似然的模型，解决剪接动力学的完整基因转录动力学，它由两组参数控制：（1）转录、剪接和降解的反应速率（2）细胞-转录状态和时间的特定潜在变量。通过 EM 迭代地推断参数。对于给定的反应速率参数估计，通过最小化每个单元与当前相轨迹的距离来将时间点分配给每个单元。转录状态是通过将可能性与轨迹上的各个片段相关联来分配的——即激活、抑制以及活跃和不活跃的稳态。d.然后通过更新反应速率的模型参数来优化整体可能性。紫色虚线将推断的（未观察到的）非活动状态与活动稳定状态联系起来。

### 2.分析结果

#### 2.1 样本spliced/unspliced的比例统计

此处展示拼接/未拼接计数的比例，通常有 10%-25% 的未拼接分子包含内含子序列。

![image-20211217174456819](C:\Users\jmzha\AppData\Roaming\Typora\typora-user-images\image-20211217174456819.png)

> 图2：spliced(mature) /unspliced(nascent) 计数比例
>
> 图注：纵坐标代表细胞类型，横坐标代表未剪接RNA（spliced，红色）和剪接mRNA（unspliced，蓝色）丰度比例
>
> `path: BMK_1_data/*_proportions`*

#### 2.2 RNA-速度

速度是基因表达空间中的向量，代表单个细胞运动的方向和速度。速度是通过模拟剪接动力学的转录动力学获得的，对于每个基因，拟合了成熟（未剪接）和成熟（剪接）mRNA 计数的稳态比率，这构成了恒定的转录状态。正速度表明基因被上调，这种情况发生在该基因未剪接 mRNA 丰度高于预期的稳定状态的细胞中。相反，负速度表明基因被下调。如图将RNA速度映射到细胞UMAP降维结果。

<img src="C:\Users\jmzha\AppData\Roaming\Typora\typora-user-images\image-20211217174046885.png" alt="image-20211217174046885" style="zoom: 67%;" />

> 图3 在单个细胞水平标明RNA-速度的方向和大小
>
> 图注：图中每个点代表单个细胞，相同颜色的点为同一类型的细胞，箭头标明了单个细胞的RNA速度的方向和大小
>
> `path: BMK_1_data/*_velocity_embedding_dy*`

<img src="C:\Users\jmzha\AppData\Roaming\Typora\typora-user-images\image-20211217174110570.png" alt="image-20211217174110570" style="zoom: 67%;" />

> 图4 单个细胞水平标明RNA-速度的方向和大小(流线图)
>
> 图注：每一种颜色代表一个细胞簇cluster，箭头代表细胞发育轨迹的方向
>
> `path: BMK_1_data/*_embedding_stream_dy*`

#### 2.3 基于RNA速度latent time

动力学模型恢复了细胞过程的潜在时间。潜在时间代表细胞的内部时钟，并近似于细胞在分化时经历的真实时间，仅基于其转录动力学。通过计算每个细胞分化需要的时间，模拟细胞发育的时间长度，以latent time图代表细胞的内部分化发育时间，近似于细胞分化时所经历需要的实时时间，从而推断细胞分化的轨迹方向，判断细胞的分化时期和细胞的根起源和终点。一般细胞分化所需要的时间也短越趋向于细胞发育的根起源，分化所需要的的时间越久越趋向于细胞分化的终点。

<img src="C:\Users\jmzha\AppData\Roaming\Typora\typora-user-images\image-20211217174127965.png" alt="image-20211217174127965" style="zoom: 67%;" />

> 图5 细胞潜在时间
>
> 图注：每个小点代表单个细胞，颜色越深代表距离开始分化时期越近，颜色越前表示距离开始分化时期越远
>
> `path: BMK_1_data/*_latent_time_dy*`

#### 2.4 驱动基因

通过RNA速度分析寻找驱动细胞RNA表达、细胞分化关键性基因，基于RNA速度分析建立了动力学模型，动态模型允许系统地将假定的驱动基因识别为具有高可能性特征的基因，高可能性选择的基因表现出明显的动态行为，而低可能性基因的表达受噪音或不存在的瞬时状态控制。这里选择RNA速度动态模型最相关TOP15驱动基因进行热图、散点图分析,揭示细胞分化的起源和终端，以及动态过程中驱动基因RNA速度变化。

<img src="C:\Users\jmzha\AppData\Roaming\Typora\typora-user-images\image-20211217174207238.png" alt="image-20211217174207238" style="zoom: 67%;" />

#### 2.5 基于RNA速度的PAGA分析

PAGA提供了轨迹推断构图的方法，其中加权边对应于两个集群之间的连接。在这里，PAGA 通过速度推断的方向性进行了扩展。

<img src="C:\Users\jmzha\AppData\Roaming\Typora\typora-user-images\image-20211217174218040.png" alt="image-20211217174218040" style="zoom: 67%;" />

> 图7 基于RNA速度的PAGA分析
>
> 图注：图中每个点代表单个细胞，相同颜色的点为同一类型的细胞，箭头代表细胞潜在的分化方向
>
> `path: BMK_1_data/*_paga_dy*`
