import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# ==============================================================================
# 1. 페이지 환경설정
# ==============================================================================
st.set_page_config(
    page_title="Ether-PC CID MS/MS Fragmentation Dashboard",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ==============================================================================
# 2. FRAGMENT_DATABASE 정의
# ==============================================================================
FRAGMENT_DATABASE = {
    "Precursor [M+Na]+": {
        "name": "Precursor Ion [M+Na]+",
        "theoretical_mz": 814.5721,
        "smiles": "C[N+](C)(C)CCOP(=O)([O-])OCC(CO[*:1])OC(=O)[*:2].[Na+]",
        "formula": "[C46H82NO7P + Na]+",
        "cleavage_mechanism": "미해리 생존 온전한 나트륨 부착 전구체 분자종",
        "reference": "Colsch et al. (Fig 1); Han & Gross (Scheme I); Al-Saad et al. (Scheme 1c)",
        "status": "Validated",
        "filter_reason": None,
        "ce_profile": {
            "30eV": {"observed_mz": 814.5746, "rel_abundance": 24.40},
            "33eV": {"observed_mz": 814.5750, "rel_abundance": 8.71},
            "37eV": {"observed_mz": 814.5774, "rel_abundance": 1.94},
            "40eV": {"observed_mz": 814.5919, "rel_abundance": 1.00},
            "50eV": {"observed_mz": 814.5225, "rel_abundance": 1.41},
        },
    },
    "[M+Na - TMA]+": {
        "name": "[M+Na - TMA]+",
        "theoretical_mz": 755.4987,
        "smiles": "O=C(OC(COP1(=O)OCCO1)CO[*:1])[*:2].[Na+]",
        "formula": "[M+Na - C3H9N]+",
        "cleavage_mechanism": "인산 음이온의 분자내 친핵 공격을 통한 5원자 환형 인산 형성 및 트리메틸아민(-59.07 Da) 탈락",
        "reference": "Colsch et al. (Scheme 6a); Han & Gross (Scheme I, Fig 2); Al-Saad et al. (Table 1 Ion B, Scheme 2)",
        "status": "Validated",
        "filter_reason": None,
        "ce_profile": {
            "30eV": {"observed_mz": 755.4996, "rel_abundance": 100.00},
            "33eV": {"observed_mz": 755.4992, "rel_abundance": 100.00},
            "37eV": {"observed_mz": 755.4993, "rel_abundance": 53.00},
            "40eV": {"observed_mz": 755.4999, "rel_abundance": 23.32},
            "50eV": {"observed_mz": 755.4677, "rel_abundance": 0.16},
        },
    },
    "[M+Na - 183]+": {
        "name": "[M+Na - 183]+",
        "theoretical_mz": 631.5051,
        "smiles": "C=C(CO[*:1])OC(=O)[*:2].[Na+]",
        "formula": "[M+Na - C5H14NO4P]+",
        "cleavage_mechanism": "sn-3 Non-sodiated Phosphocholine(-183.07 Da) 중성 탈락 (소듐 유지 에놀 에테르/에스터)",
        "reference": "Godzien et al. (Table 1); Al-Saad et al. (Table 1 Ion C, Scheme 2); Han & Gross (Scheme I Pathway 4)",
        "status": "Validated",
        "filter_reason": None,
        "ce_profile": {
            "30eV": {"observed_mz": 631.5068, "rel_abundance": 36.22},
            "33eV": {"observed_mz": 631.5070, "rel_abundance": 74.76},
            "37eV": {"observed_mz": 631.5065, "rel_abundance": 100.00},
            "40eV": {"observed_mz": 631.5063, "rel_abundance": 100.00},
            "50eV": {"observed_mz": 631.5108, "rel_abundance": 6.55},
        },
    },
    "[M+Na - 205]+": {
        "name": "[M+Na - 205]+",
        "theoretical_mz": 609.5242,
        "smiles": "[*:1]OCC1C[O+]=C([*:2])O1",
        "formula": "[M+Na - C5H14NO4PNa]+",
        "cleavage_mechanism": "sn-3 Sodiated Phosphocholine(-205.05 Da) 중성 탈락 및 sn-2 아실 카보닐의 백본 분자내 고리화(1,3-dioxolan-2-ylium 양이온 형성)",
        "reference": "Han & Gross (Fig 2b/c, Scheme I Pathway 3); Al-Saad et al. (Table 1 Ion D, Scheme 5); Colsch et al. (Table 1a)",
        "status": "Validated",
        "filter_reason": None,
        "ce_profile": {
            "30eV": {"observed_mz": 609.5255, "rel_abundance": 7.16},
            "33eV": {"observed_mz": 609.5262, "rel_abundance": 12.84},
            "37eV": {"observed_mz": 609.5261, "rel_abundance": 13.94},
            "40eV": {"observed_mz": 609.5264, "rel_abundance": 8.39},
            "50eV": {"observed_mz": 609.5359, "rel_abundance": 0.37},
        },
    },
    "Protonated Phosphocholine": {
        "name": "Protonated Phosphocholine",
        "theoretical_mz": 184.0733,
        "smiles": "C[N+](C)(C)CCOP(=O)(O)O",
        "formula": "C5H15NO4P+",
        "cleavage_mechanism": "sn-3 Phosphodiester 결합 해리를 통한 포스포콜린 극성 머리그룹 양이온 형성",
        "reference": "Colsch et al. (Table 1a, Fig 3d); Godzien et al. (Table 1); Han & Gross (Fig 2); Al-Saad et al. (Fig 2, Scheme 8)",
        "status": "Validated",
        "filter_reason": None,
        "ce_profile": {
            "30eV": {"observed_mz": 184.0744, "rel_abundance": 0.64},
            "33eV": {"observed_mz": 184.0742, "rel_abundance": 0.67},
            "37eV": {"observed_mz": 184.0738, "rel_abundance": 2.42},
            "40eV": {"observed_mz": 184.0733, "rel_abundance": 2.53},
            "50eV": {"observed_mz": 184.0744, "rel_abundance": 6.34},
        },
    },
    "Sodiated Cyclophosphane": {
        "name": "Sodiated Cyclophosphane",
        "theoretical_mz": 146.9817,
        "smiles": "O=P1(O)OCCO1.[Na+]",
        "formula": "C2H4O4PNa+",
        "cleavage_mechanism": "헤드그룹 유래 5원자 환형 고리형 인산(1,3,2-dioxaphospholane 2-oxide) 나트륨 착이온",
        "reference": "Colsch et al. (Scheme 6c/d); Godzien et al. (Table 1); Han & Gross (Scheme I); Al-Saad et al. (Table 1 Ion L, Scheme 2)",
        "status": "Validated",
        "filter_reason": None,
        "ce_profile": {
            "30eV": {"observed_mz": 146.9829, "rel_abundance": 3.94},
            "33eV": {"observed_mz": 146.9819, "rel_abundance": 13.63},
            "37eV": {"observed_mz": 146.9817, "rel_abundance": 32.15},
            "40eV": {"observed_mz": 146.9823, "rel_abundance": 59.85},
            "50eV": {"observed_mz": 146.9822, "rel_abundance": 100.00},
        },
    },
    "Protonated Ethylene Phosphate": {
        "name": "Protonated Ethylene Phosphate",
        "theoretical_mz": 125.0009,
        "smiles": "O=P1(O)OCCO1",
        "formula": "C2H6O4P+",
        "cleavage_mechanism": "포스포콜린 머리그룹에서 4차 암모늄 탈락 후 잔존한 5원자 환형 인산 에스테르 양이온",
        "reference": "Colsch et al. (Table 1a, Scheme 5b); Godzien et al. (Table 1)",
        "status": "Excluded",
        "filter_reason": "Low Intensity (<=1% at all CE)",
        "ce_profile": {
            "30eV": {"observed_mz": 125.1370, "rel_abundance": 0.01},
            "33eV": {"observed_mz": 125.1348, "rel_abundance": 0.07},
            "37eV": {"observed_mz": 125.1247, "rel_abundance": 0.04},
            "40eV": {"observed_mz": 125.1343, "rel_abundance": 0.32},
            "50eV": {"observed_mz": 125.0021, "rel_abundance": 1.00},
        },
    },
    "Protonated Choline": {
        "name": "Protonated Choline",
        "theoretical_mz": 104.1070,
        "smiles": "C[N+](C)(C)CCO",
        "formula": "C5H14NO+",
        "cleavage_mechanism": "5가 인 중간체 형성을 거친 분자내 수소 전이 기반의 콜린 양이온 형성",
        "reference": "Colsch et al. (Table 1a, Scheme 6b)",
        "status": "Validated",
        "filter_reason": None,
        "ce_profile": {
            "30eV": {"observed_mz": 104.1084, "rel_abundance": 0.20},
            "33eV": {"observed_mz": 104.1077, "rel_abundance": 0.37},
            "37eV": {"observed_mz": 104.1074, "rel_abundance": 1.58},
            "40eV": {"observed_mz": 104.1071, "rel_abundance": 4.98},
            "50eV": {"observed_mz": 104.1070, "rel_abundance": 14.67},
        },
    },
    "Trimethyl vinyl ammonium": {
        "name": "Trimethyl vinyl ammonium",
        "theoretical_mz": 86.0964,
        "smiles": "C=C[N+](C)(C)C",
        "formula": "C5H12N+",
        "cleavage_mechanism": "헤드그룹 인산 에스테르 C-O 결합 절단 및 비닐기 전이를 수반한 4차 암모늄 형성",
        "reference": "Colsch et al. (Scheme 5a); Al-Saad et al. (Table 1 Ion M, Scheme 7)",
        "status": "Validated",
        "filter_reason": None,
        "ce_profile": {
            "30eV": {"observed_mz": 86.0972, "rel_abundance": 3.20},
            "33eV": {"observed_mz": 86.0967, "rel_abundance": 5.63},
            "37eV": {"observed_mz": 86.0970, "rel_abundance": 8.70},
            "40eV": {"observed_mz": 86.0969, "rel_abundance": 10.47},
            "50eV": {"observed_mz": 86.0971, "rel_abundance": 29.23},
        },
    },
    "Odd-electron Nitrogen ion": {
        "name": "Odd-electron Nitrogen ion",
        "theoretical_mz": 71.0730,
        "smiles": "C[N+]1(C)C[CH]C1",
        "formula": "C4H9N+•",
        "cleavage_mechanism": "헤드그룹 4차 암모늄 부위 C-N 라디칼 절단 (dehydrogenated 1,1-dimethylazetidinium)",
        "reference": "Colsch et al. (Scheme 5a)",
        "status": "Excluded",
        "filter_reason": "LPC literature origin (Incompatible with Ether-PC ESI)",
        "ce_profile": {
            "30eV": {"observed_mz": 71.0860, "rel_abundance": 0.30},
            "33eV": {"observed_mz": 71.0857, "rel_abundance": 1.24},
            "37eV": {"observed_mz": 71.0857, "rel_abundance": 3.58},
            "40eV": {"observed_mz": 71.0855, "rel_abundance": 4.58},
            "50eV": {"observed_mz": 71.0862, "rel_abundance": 10.42},
        },
    },
}

# ==============================================================================
# 3. RDKit 2D 렌더링 헬퍼 함수
# ==============================================================================
def get_mol_image(smiles_str: str, size=(300, 220)) -> Image.Image:
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "*":
            if atom.GetAtomMapNum() == 1:
                atom.SetProp("atomLabel", "R1")
            elif atom.GetAtomMapNum() == 2:
                atom.SetProp("atomLabel", "R2")
            else:
                atom.SetProp("atomLabel", "R")
    return Draw.MolToImage(mol, size=size)

# ==============================================================================
# 4. 사이드바 인터페이스
# ==============================================================================
with st.sidebar:
    st.header("⚙️ 분석 파라미터")
    
    tolerance = st.slider(
        "허용 오차 (Tolerance, ±Da)",
        min_value=0.0000,
        max_value=0.1000,
        value=0.0100,
        step=0.0005,
        format="%.4f",
        help="슬라이더를 0으로 줄이면 오차가 큰 피크가 실시간으로 표와 그리드에서 제외됩니다."
    )
    
    st.markdown("---")
    validation_filter = st.toggle(
        "연구자 검증 필터 적용 (Researcher Validation)",
        value=True,
        help="ON: 제외 사유가 있는 2종(m/z 125, m/z 71)을 제외하고 8종만 표시합니다."
    )
    
    st.info(
        "💡 **필터 동작 기준**\n\n"
        f"- 현재 허용 오차: **±{tolerance:.4f} Da**\n"
        f"- 필터 모드: **{'검증 완료 (8종)' if validation_filter else '전체 표시 (10종)'}**"
    )

# ==============================================================================
# 5. 실시간 데이터 필터링 & DataFrame 생성
# ==============================================================================
filtered_items = {}
table_rows = []

for key, item in FRAGMENT_DATABASE.items():
    # 1) 검증 필터 적용 여부
    if validation_filter and item["status"] == "Excluded":
        continue
    
    # 2) Tolerance 실시간 필터 적용 (30eV 기준 오차)
    obs_30 = item["ce_profile"]["30eV"]["observed_mz"]
    delta_da = obs_30 - item["theoretical_mz"]
    
    if abs(delta_da) <= tolerance:
        filtered_items[key] = item
        table_rows.append({
            "Fragment Name": item["name"],
            "Formula": item["formula"],
            "Theoretical m/z": f"{item['theoretical_mz']:.4f}",
            "Observed m/z (30eV)": f"{obs_30:.4f}",
            "Delta (Da)": f"{delta_da:+.4f}",
            "Status": "✅ " + item["status"] if item["status"] == "Validated" else "⚠️ Excluded",
            "Reason / Reference": item["filter_reason"] if item["filter_reason"] else item["reference"].split(";")[0],
        })

df_summary = pd.DataFrame(table_rows)

# ==============================================================================
# 6. 메인 화면 구성
# ==============================================================================
st.title("🧪 Ether-PC CID MS/MS Fragmentation Dashboard")
st.caption("Collision-Induced Dissociation Structural Annotation Platform for Ether-linked Phosphatidylcholines")

# 타겟 지질 헤더 요약 정보
header_col1, header_col2, header_col3 = st.columns([4, 4, 3])
with header_col1:
    st.markdown("**Target Lipid**: Ether-linked PC (Plasmanyl/Plasmenyl)")
with header_col2:
    st.markdown("**Precursor Adduct**: `[M+Na]+` (m/z 814.5721)")
with header_col3:
    st.markdown(f"**Matched Fragments**: `{len(filtered_items)} / 10 Active`")

st.markdown("---")

tab1, tab2 = st.tabs(["🔬 Fragment 매칭 & RDKit 2D 화학 구조", "📈 충돌 에너지(CE) 거동 분석"])

# ------------------------------------------------------------------------------
# 탭 1: Fragment 매칭 요약표 & RDKit 2D 구조 그리드
# ------------------------------------------------------------------------------
with tab1:
    st.subheader("📋 실시간 Fragment 매칭 결과 요약")
    
    if not df_summary.empty:
        st.dataframe(df_summary, use_container_width=True, hide_index=True)
    else:
        st.warning(f"현재 설정된 허용 오차(±{tolerance:.4f} Da) 내에 매칭되는 Fragment가 없습니다. 슬라이더 값을 올려주세요.")
    
    st.markdown("---")
    st.subheader("🧩 RDKit 2D 화학 구조 분해 카드 그리드")
    
    if filtered_items:
        cols = st.columns(3)
        for idx, (k, item) in enumerate(filtered_items.items()):
            col = cols[idx % 3]
            with col:
                with st.container(border=True):
                    # 헤더 및 m/z
                    st.markdown(f"**{item['name']}**")
                    st.caption(f"이론 m/z: **{item['theoretical_mz']:.4f}** | 식: `{item['formula']}`")
                    
                    # m/z 609.5242 검증 배지 표출
                    if "609" in k:
                        st.success("🏷️ [Validated Structure: ESI-specific]")
                    
                    # RDKit 2D 렌더링
                    img = get_mol_image(item["smiles"])
                    if img is not None:
                        st.image(img, use_container_width=True)
                    else:
                        st.warning("화학 구조 렌더링 불가")
                    
                    # 기작 및 출처
                    st.markdown(f"**개열 기작**: {item['cleavage_mechanism']}")
                    st.caption(f"**출처**: {item['reference']}")
                    
                    if item["status"] == "Excluded":
                        st.error(f"제외 사유: {item['filter_reason']}")
    else:
        st.info("표시할 단편 구조 카드가 없습니다.")

# ------------------------------------------------------------------------------
# 탭 2: 충돌 에너지(CE: 30~50 eV) 거동 분석 인터랙티브 차트
# ------------------------------------------------------------------------------
with tab2:
    st.subheader("📈 Collision Energy (30 ~ 50 eV) Breakdown Curves")
    st.markdown("충돌 에너지 증가에 따른 단편 이온들의 상대강도(%) 변화를 추적합니다.")
    
    available_keys = list(filtered_items.keys())
    
    if available_keys:
        default_targets = [
            "Precursor [M+Na]+",
            "[M+Na - TMA]+",
            "[M+Na - 183]+",
            "[M+Na - 205]+",
            "Sodiated Cyclophosphane"
        ]
        default_selected = [t for t in default_targets if t in available_keys]
        
        selected_fragments = st.multiselect(
            "차트에 표시할 Fragment 선택:",
            options=available_keys,
            default=default_selected if default_selected else available_keys[:3]
        )
        
        if selected_fragments:
            ce_labels = ["30eV", "33eV", "37eV", "40eV", "50eV"]
            fig = go.Figure()
            
            for frag_key in selected_fragments:
                frag_data = filtered_items[frag_key]
                y_values = [frag_data["ce_profile"][ce]["rel_abundance"] for ce in ce_labels]
                
                fig.add_trace(go.Scatter(
                    x=ce_labels,
                    y=y_values,
                    mode="lines+markers",
                    name=f"{frag_data['name']} (m/z {frag_data['theoretical_mz']:.1f})",
                    line=dict(width=2.5),
                    marker=dict(size=7)
                ))
            
            fig.update_layout(
                title="Relative Abundance (%) vs Collision Energy (eV)",
                xaxis_title="Collision Energy (CE)",
                yaxis_title="Relative Abundance (% of Base Peak)",
                template="plotly_white",
                hovermode="x unified",
                height=520,
                legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("그래프에 표시할 Fragment를 1개 이상 선택하세요.")
    else:
        st.warning("선택 가능한 활성 Fragment가 없습니다.")

# ==============================================================================
# 7. 결과 내보내기 (CSV 다운로드 버튼)
# ==============================================================================
st.markdown("---")
col_exp1, col_exp2 = st.columns([8, 2])

with col_exp1:
    st.caption("현재 화면에 필터링된 최종 매칭 결과를 UTF-8-BOM CSV 파일로 내보냅니다.")

with col_exp2:
    if not df_summary.empty:
        csv_data = df_summary.to_csv(index=False, encoding="utf-8-sig")
        st.download_button(
            label="📥 CSV 결과 다운로드",
            data=csv_data,
            file_name=f"EtherPC_Annotation_tol_{tolerance:.4f}.csv",
            mime="text/csv",
            use_container_width=True
        )
    else:
        st.button("📥 CSV 결과 다운로드", disabled=True, use_container_width=True)