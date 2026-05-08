import {
  BarChart,
  Callout,
  Card,
  CardBody,
  CardHeader,
  Divider,
  Grid,
  H1,
  H2,
  H3,
  Pill,
  Row,
  Stack,
  Stat,
  Table,
  Text,
  useHostTheme,
} from "cursor/canvas";

const runRows = [
  {
    runId: "R001",
    rankingCount: 12,
    selectedCount: 2,
    mrnaLength: 193,
    gcPercent: 26.42,
    mfe: -42.6,
    status: "阶段验收通过",
  },
  {
    runId: "R002",
    rankingCount: 18,
    selectedCount: 3,
    mrnaLength: 229,
    gcPercent: 35.81,
    mfe: -72.6,
    status: "阶段验收通过",
  },
  {
    runId: "R003",
    rankingCount: 18,
    selectedCount: 3,
    mrnaLength: 229,
    gcPercent: 35.81,
    mfe: -72.6,
    status: "阶段验收通过",
  },
  {
    runId: "R_public_001",
    rankingCount: 198,
    selectedCount: 10,
    mrnaLength: 481,
    gcPercent: 51.35,
    mfe: -274.5,
    status: "阶段验收通过",
  },
];

const progressRows = [
  ["输入契约", "已完成", "BioDriver 候选新抗原、HLA 分型、meta 已接入并可校验"],
  ["主流程", "已完成", "输入校验 -> 评分排序 -> 肽段筛选 -> mRNA 构建 -> QC 报告 -> SimHub 交付"],
  ["真实后端", "已完成", "NetMHCpan、NetMHCIIpan、DeepImmuno、PRIME、Repitope 已有真实来源口径"],
  ["结构交付", "已完成", "4 个实例均生成 PANDORA 真实 peptide-MHC 结构 complex.pdb"],
  ["mRNA 设计", "已完成", "LinearDesign 真实命令行接入，4 个实例均完成密码子优化"],
  ["稳定性质控", "已完成", "RNAfold、RNAeval、RNAplfold 均正常输出可复检指标"],
  ["SimHub 回传", "等待下游", "输入包已交付，仍待 MD 轨迹、能量、RMSD 和 QC flags 回传"],
  ["临床终版", "等待上游", "R002/R003/R_public_001 仍需替换为 BioDriver 患者真实 HLA-II 分型后重跑"],
];

const deliverables = [
  "4 个阶段验收实例：R001、R002、R003、R_public_001",
  "每个实例均有 mRNA FASTA、设计元数据、稳定性质控 JSON、报告和自证材料",
  "每个实例均有 SimHub 交付包，包含 PANDORA 真实 complex.pdb",
  "关键交付文件未发现 proxy、fallback、coarse 残留",
];

const modelStageRows = [
  ["候选输入", "BioDriver", "提供患者突变肽、VAF、HLA 分型", "上游数据源，不是预测模型"],
  ["MHC-I 结合", "NetMHCpan", "判断肽段是否容易被患者 I 类 HLA 呈递", "已接入真实后端"],
  ["MHC-I 增强", "BigMHC", "作为第二路交叉验证，减少单模型偏差", "可选增强，当前默认关闭"],
  ["MHC-II 结合", "NetMHCIIpan", "覆盖 CD4 T 细胞相关呈递信号", "已接入，终版等待真实 HLA-II"],
  ["免疫原性", "DeepImmuno / PRIME / Repitope", "判断肽段更可能激活免疫反应", "已用真实来源口径"],
  ["结构建模", "PANDORA + MODELLER", "生成 peptide-MHC 真实结构给 SimHub", "4 个实例已完成"],
  ["结构复核", "ColabFold / AFM", "对重点候选做独立结构复核", "可选，成本更高（本轮 AFM 置信度提升）"],
  ["mRNA 设计", "LinearDesign", "进行密码子优化并生成 mRNA 序列", "当前主力方案"],
  ["mRNA 稳定性", "RNAfold / RNAeval / RNAplfold", "检查二级结构和稳定性风险", "4 个实例均通过"],
  ["下游验证", "SimHub MD", "用分子动力学看结构稳定性", "等待回传证据"],
];

const mrnaStrategyRows = [
  ["A 稳妥版", "沿用现有 Top 肽 + LinearDesign", "最低", "最快形成交付，适合当前验收版", "已完成"],
  ["B 免疫优先版", "提高 DeepImmuno / PRIME / Repitope 权重", "低-中", "更重视潜在免疫反应，可能牺牲部分稳定性", "建议新增"],
  ["C 结构优先版", "PANDORA 通过后再用 ColabFold/AFM 复核 Top 肽", "中-高", "减少结构不稳定候选，适合进 SimHub 前筛选", "建议新增"],
  ["D 稳定性优先版", "同一肽组合生成多条 mRNA，按 RNAfold/RNAplfold 选优", "中", "提高表达与稳定性把握，适合最终候选收敛", "建议新增"],
  ["E 成本优先版", "只保留多模型一致高分候选，减少 SimHub 数量", "最低", "省钱但探索空间较小，适合预算受限场景", "建议新增"],
];

function StatusLine({ title, body }: { title: string; body: string }) {
  const theme = useHostTheme();

  return (
    <div
      style={{
        borderLeft: `3px solid ${theme.accent.primary}`,
        paddingLeft: 12,
        paddingTop: 2,
        paddingBottom: 2,
      }}
    >
      <Text weight="semibold">{title}</Text>
      <Text tone="secondary" size="small">
        {body}
      </Text>
    </div>
  );
}

export default function ImmunoGenExecutiveReport() {
  const theme = useHostTheme();
  const totalRanking = runRows.reduce((sum, item) => sum + item.rankingCount, 0);
  const totalSelected = runRows.reduce((sum, item) => sum + item.selectedCount, 0);

  return (
    <Stack gap={22} style={{ maxWidth: 1180 }}>
      <Stack gap={10}>
        <Row justify="space-between" align="start" gap={16} wrap>
          <Stack gap={6} style={{ maxWidth: 760 }}>
            <H1>ImmunoGen 项目报告</H1>
            <Text tone="secondary">
              目标：把 BioDriver 提供的患者新抗原和 HLA 分型，转成可交付的个性化多价 mRNA 疫苗设计，并准备 SimHub 分子动力学验证包。
            </Text>
          </Stack>
          <Row gap={8} wrap>
            <Pill tone="success" active>
              阶段验收版
            </Pill>
            <Pill tone="info">最近核对：2026-04-30</Pill>
          </Row>
        </Row>

        <Callout tone="success" title="一句话结论">
          项目主链路已经从“能跑”推进到“可验收、可复检、可交付”：4 个实例已完成真实后端严格重跑。下一步建议用多模型组合生成 3-5 条 mRNA 候选，再按效果、稳定性、验证成本做性价比选择。
        </Callout>
      </Stack>

      <Grid columns={4} gap={14}>
        <Stat value="4/4" label="实例阶段验收通过" tone="success" />
        <Stat value={String(totalRanking)} label="候选排序记录" tone="info" />
        <Stat value={String(totalSelected)} label="入选疫苗肽段" tone="success" />
        <Stat value="3-5" label="建议并行 mRNA 候选" tone="warning" />
      </Grid>

      <Grid columns="1.15fr 0.85fr" gap={18} align="stretch">
        <Stack gap={12}>
          <H2>老板关心的项目状态</H2>
          <Grid columns={2} gap={14}>
            <StatusLine
              title="能不能交付？"
              body="可以交付阶段验收版和下游技术检查版；4 个 run 均有报告、自证、mRNA 设计和 SimHub 输入包。"
            />
            <StatusLine
              title="是不是真实工具？"
              body="核心后端已切换到真实来源，并有防回退检查；交付包未发现 proxy、fallback、coarse 残留。"
            />
            <StatusLine
              title="还差什么？"
              body="差临床终版 HLA-II 数据和 SimHub MD 回传结果，这两项分别依赖 BioDriver 和 SimHub。"
            />
            <StatusLine
              title="风险在哪里？"
              body="计算设计不能替代湿实验；MHC-II 终版数据到位后，需要自动重跑并更新全部验收材料。"
            />
            <StatusLine
              title="为什么要多模型？"
              body="单一模型容易偏向某一类指标；多模型一致高分的候选，更适合进入昂贵的 SimHub 和实验验证。"
            />
            <StatusLine
              title="怎么选性价比？"
              body="先低成本批量筛选，再把少数高分候选送结构和 MD，最后按通过率、成本和交付速度选 1-2 条主线。"
            />
          </Grid>
        </Stack>

        <Card>
          <CardHeader trailing={<Pill tone="warning" size="sm">待外部输入</Pill>}>当前阻塞</CardHeader>
          <CardBody>
            <Stack gap={12}>
              <Text>
                <Text weight="semibold" as="span">BioDriver：</Text>
                需要提供 R002、R003、R_public_001 的患者真实 HLA-II 临床终版分型。
              </Text>
              <Divider />
              <Text>
                <Text weight="semibold" as="span">SimHub：</Text>
                需要回传 MD 轨迹、能量、RMSD、QC flags 和 summary，完成稳定性证据归档。
              </Text>
            </Stack>
          </CardBody>
        </Card>
      </Grid>

      <H2>每个阶段使用的模型</H2>
      <Text tone="secondary">
        管理层理解口径：这些模型不是都直接生成 mRNA，而是分工合作。前面模型负责挑“值得做”的肽，中间模型负责看结构是否靠谱，后面模型负责把肽组合成 mRNA 并检查稳定性。
      </Text>
      <Table
        headers={["阶段", "模型/工具", "老板能理解的作用", "当前状态"]}
        rows={modelStageRows}
        rowTone={modelStageRows.map((row) => {
          if (row[3] === "已接入真实后端" || row[3] === "4 个实例已完成" || row[3] === "当前主力方案" || row[3] === "4 个实例均通过") return "success";
          if (row[3].includes("等待") || row[3].includes("可选")) return "warning";
          return "neutral";
        })}
        striped
      />

      <Grid columns="0.9fr 1.1fr" gap={18} align="stretch">
        <Card>
          <CardHeader trailing={<Pill tone="info" size="sm">建议策略</Pill>}>多 mRNA 生成思路</CardHeader>
          <CardBody>
            <Stack gap={12}>
              <Text>
                建议不要只押一条 mRNA。更稳妥的方式是让同一批候选肽经过不同模型权重和不同优化目标，形成 3-5 条候选 mRNA。
              </Text>
              <Divider />
              <Text tone="secondary">
                第一轮用低成本模型快速筛掉明显不合适的候选；第二轮只把少量高分候选送结构复核和 SimHub；第三轮结合成本、时间和稳定性，选 1-2 条进入更深验证。
              </Text>
            </Stack>
          </CardBody>
        </Card>

        <Stack gap={12}>
          <H2>候选 mRNA 性价比方案</H2>
          <Table
            headers={["方案", "生成方式", "成本", "适用场景", "状态"]}
            rows={mrnaStrategyRows}
            rowTone={mrnaStrategyRows.map((row) => (row[4] === "已完成" ? "success" : "info"))}
            striped
          />
        </Stack>
      </Grid>

      <H2>4 个实例的关键数据</H2>
      <Grid columns="1fr 1fr" gap={18}>
        <Card>
          <CardHeader>候选量与入选量</CardHeader>
          <CardBody>
            <BarChart
              categories={runRows.map((item) => item.runId)}
              series={[
                { name: "候选排序记录", data: runRows.map((item) => item.rankingCount), tone: "info" },
                { name: "入选肽段", data: runRows.map((item) => item.selectedCount), tone: "success" },
              ]}
              height={260}
            />
          </CardBody>
        </Card>

        <Card>
          <CardHeader>mRNA 长度与 GC 含量</CardHeader>
          <CardBody>
            <BarChart
              categories={runRows.map((item) => item.runId)}
              series={[
                { name: "mRNA 长度", data: runRows.map((item) => item.mrnaLength), tone: "neutral" },
                { name: "GC 含量", data: runRows.map((item) => item.gcPercent), tone: "warning" },
              ]}
              height={260}
            />
          </CardBody>
        </Card>
      </Grid>

      <Table
        headers={["实例", "排序记录", "入选肽段", "mRNA 长度", "GC 含量", "RNAfold MFE", "状态"]}
        rows={runRows.map((item) => [
          item.runId,
          item.rankingCount,
          item.selectedCount,
          item.mrnaLength,
          `${item.gcPercent}%`,
          item.mfe,
          item.status,
        ])}
        columnAlign={["left", "right", "right", "right", "right", "right", "left"]}
        rowTone={runRows.map(() => "success")}
        striped
      />

      <Grid columns="0.95fr 1.05fr" gap={18}>
        <Stack gap={12}>
          <H2>已形成的交付物</H2>
          <Card>
            <CardHeader>可给下游检查的内容</CardHeader>
            <CardBody>
              <Stack gap={10}>
                {deliverables.map((item) => (
                  <Row key={item} gap={10} align="start">
                    <span style={{ color: theme.accent.primary, lineHeight: "20px" }}>•</span>
                    <Text>{item}</Text>
                  </Row>
                ))}
              </Stack>
            </CardBody>
          </Card>
        </Stack>

        <Stack gap={12}>
          <H2>模块完成情况</H2>
          <Table
            headers={["模块", "状态", "说明"]}
            rows={progressRows}
            rowTone={progressRows.map((row) => {
              if (row[1] === "已完成") return "success";
              if (row[1] === "等待下游" || row[1] === "等待上游") return "warning";
              return "neutral";
            })}
            striped
          />
        </Stack>
      </Grid>

      <Divider />

      <Grid columns={3} gap={16}>
        <Card>
          <CardHeader>下一步 1</CardHeader>
          <CardBody>
            <H3>接收 SimHub 回传</H3>
            <Text tone="secondary">
              收到轨迹、能量、RMSD 和 QC flags 后，归档为 evidence，并更新 SELF_CHECK 与 SIMHUB_EVIDENCE。
            </Text>
          </CardBody>
        </Card>
        <Card>
          <CardHeader>下一步 2</CardHeader>
          <CardBody>
            <H3>替换临床 HLA-II</H3>
            <Text tone="secondary">
              BioDriver 数据到位后，对 R002、R003、R_public_001 执行 MHC-II 到 SimHub 的全链路重跑。
            </Text>
          </CardBody>
        </Card>
        <Card>
          <CardHeader>下一步 3</CardHeader>
          <CardBody>
            <H3>清理提交边界</H3>
            <Text tone="secondary">
              提交前排除本机缓存、外部工具目录和大型数据库，避免把 .autodl、external_refs、下载缓存误提交。
            </Text>
          </CardBody>
        </Card>
      </Grid>

      <Text tone="tertiary" size="small">
        口径说明：本报告根据 README、docs/TODO.md、4 个 results/REPORT.md 与 qc_metrics.json 汇总；计算验证用于阶段筛选和下游模拟准备，不等同于临床疗效结论。
      </Text>
    </Stack>
  );
}
