import pandas as pd
import networkx as nx
import random

# 读取 Excel 文件并构建无向图
file_path = r"D:\建信金科数据集\参数\汇总小编号量化参数_更新.xlsx"
immunity_file = r"D:\建信金科数据集\参数\直接免疫率.xlsx"
output_file = r"D:\建信金科数据集\参数\风险源\感染10%终止子网络.xlsx"

data = pd.read_excel(file_path)
immunity_data = pd.read_excel(immunity_file)

# 提取免疫率数据
immunity_rates = dict(zip(immunity_data['NodeID'], immunity_data['ImmunityRate']))

# 提取 gamma 数据
gamma_values = dict(zip(data['GuarantorID'], data['immunisation_rate']))  # 以 GuarantorID 为键
default_gamma = 0.2  # 默认 gamma 值

# 构建图
G = nx.Graph()
edges = list(zip(data['GuarantorID'], data['GuaranteeID']))
G.add_edges_from(edges)

# 参数设置
beta = 0.9  # 感染率
gamma_base = 0.2  # 基础恢复率
alpha = 0.1  # 免疫丧失率
max_infected_limit = 1000  # 最大感染节点数限制

# 初始化节点状态 (S, I, R)
states = {node: 'S' for node in G.nodes}
initial_infected = [61]  # 替换为初始感染点编号
for node in initial_infected:
    states[node] = 'I'  # 初始感染点设为感染状态


# 找到初始传染源所在的子网络
for component in nx.connected_components(G):
    if 61 in component:
        subnetwork_nodes = set(component)
        break  # 找到就停止
subnetwork_size = len(subnetwork_nodes)
threshold = min(max(3, int(0.01* subnetwork_size)),100) # 30% 的节点数，最小为 3

print(f"初始传染源子网络大小: {subnetwork_size}, 终止感染阈值: {threshold}")

# 模拟传播过程
def simulate_sirs_until_stable_or_limit(graph, states, beta, gamma_base, alpha, max_infected_limit):
    time_series = {}
    all_infected_nodes = set()
    previous_infected_nodes = set()
    t = 0

    time_series[t] = list(states.items())  # 初始状态

    while True:
        new_states = states.copy()
        current_infected_nodes = set()

        for node in graph.nodes:
            if states[node] == 'S':
                # 判断是否会直接免疫
                immunity_rate = immunity_rates.get(node, 0)
                if random.random() < immunity_rate:
                    continue

                # 易感者可能被感染
                infected_neighbors = [n for n in graph.neighbors(node) if states[n] == 'I']

                # 计算 beta_i
                beta_i = 0
                for neighbor in infected_neighbors:
                    edge_data = data[
                        ((data['GuarantorID'] == neighbor) & (data['GuaranteeID'] == node)) |
                        ((data['GuarantorID'] == node) & (data['GuaranteeID'] == neighbor))
                    ]
                    if not edge_data.empty:
                        amount_num = edge_data['AmountNum'].iloc[0]
                        risk_probability = edge_data['Risk_Probability'].iloc[0]
                        base_beta_i = (amount_num * 0.5 + risk_probability * 0.5) * beta

                        # 传播方向判断
                        if ((data['GuarantorID'] == node) & (data['GuaranteeID'] == neighbor)).any():
                            beta_i += base_beta_i * 0.2  # GuaranteeID 传给 GuarantorID，降低传播率
                        else:
                            beta_i += base_beta_i  # GuarantorID 传给 GuaranteeID，保持原传播率

                beta_i = min(beta_i, 1)
                if infected_neighbors and random.random() < 1 - (1 - beta_i) ** len(infected_neighbors):
                    new_states[node] = 'I'
                    all_infected_nodes.add(node)
                    current_infected_nodes.add(node)
            elif states[node] == 'I':
                # 计算节点特定的 gamma
                node_gamma = gamma_base * gamma_values.get(node, default_gamma)
                current_infected_nodes.add(node)
                if random.random() < node_gamma:
                    new_states[node] = 'R'
            elif states[node] == 'R':
                if random.random() < alpha:
                    new_states[node] = 'S'

        # 终止条件
        if new_states == states and current_infected_nodes == previous_infected_nodes:
            break
        if len(all_infected_nodes) > max_infected_limit:
            print("Simulation stopped: Infected nodes exceeded the limit.")
            break

        previous_infected_nodes = current_infected_nodes
        states = new_states
        t += 1
        time_series[t] = list(states.items())

    return time_series, all_infected_nodes

# 解析 Excel 并转换为 Markdown 格式
def excel_to_markdown(df):
    df = df.iloc[:, :6]  # 仅取前6列
    markdown_str = "| " + " | ".join(df.columns) + " |\n"
    markdown_str += "| " + " | ".join(["-" * len(col) for col in df.columns]) + " |\n"
    for _, row in df.iterrows():
        markdown_str += "| " + " | ".join(map(str, row)) + " |\n"
    return markdown_str

# 运行模拟
infection_result, all_infected_nodes = simulate_sirs_until_stable_or_limit(
    G, states, beta, gamma_base, alpha, max_infected_limit
)
# 记录终止时刻的感染节点
infected_subnetwork_nodes = []

for time, state_list in infection_result.items():
    infected_nodes = [node for node, state in state_list if state == 'I']
    print(f"Time {time}: Infected nodes = {infected_nodes}")

    # 终止条件：感染的节点数超过阈值
    if len(infected_nodes) > threshold:
        infected_subnetwork_nodes = infected_nodes
        print(f"Simulation stopped: Infected nodes reached {len(infected_nodes)} (>{threshold}).")
        break

infected_subnetwork = data[
    (data['GuarantorID'].isin(infected_nodes)) &
    (data['GuaranteeID'].isin(infected_nodes))  # 确保两端都被感染
]

# 保存到 Excel
infected_subnetwork.to_excel(output_file, index=False)

print(f"感染10%终止子网络已保存至: {output_file}")

# 生成 Markdown 表格
markdown_table = excel_to_markdown(infected_subnetwork)

# 输出 Markdown 格式
print("\n**请分析以下感染子网络数据：**\n")
print("请判断担保网络中的风险传染源节点，答案仅输出风险源的名称即可,所有已经感染的节点如下：\n"+markdown_table)
