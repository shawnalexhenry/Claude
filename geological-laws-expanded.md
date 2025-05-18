# Expanded Geological and Geophysical Laws for Terrain Simulation

## 1. Impact Processes (Asteroid/Comet/Meteorite) - Expanded

### 1.1 Impactor Characteristics by Type

#### Asteroids
- **Density Ranges**:
  - C-type (carbonaceous): 1300-1700 kg/m³
  - S-type (silicaceous): 2200-2900 kg/m³
  - M-type (metallic): 4000-5500 kg/m³
  - Average NEO density: 2600 kg/m³
- **Typical Velocities**: 12-20 km/s (Earth impacts)
- **Shapes**: Irregular, often approximated as triaxial ellipsoids
  - Elongation ratio (a:b:c) typically 1.2:1.1:1
  - Larger asteroids (>100 km) trend toward more spherical shapes
- **Porosity**: 
  - C-type: 30-40%
  - S-type: 10-30%
  - M-type: 5-20%
- **Strength**:
  - Tensile strength: 0.01-50 MPa (highly variable)
  - Compressive strength: 1-100 MPa

#### Comets
- **Density**: 400-1200 kg/m³ (nucleus)
- **Composition**: 
  - 80-90% ice (H₂O, CO₂, CO, CH₄)
  - 10-20% refractory materials (silicates, organics)
- **Typical Velocities**: 30-70 km/s (Earth impacts)
- **Nucleus Shapes**: Highly irregular, often bilobal
  - Aspect ratios commonly 2:1:1
- **Porosity**: 40-80%
- **Strength**:
  - Extremely low cohesive strength: 0.001-0.01 MPa
  - Behave as strengthless rubble piles at scales >100 m

#### Meteorites
- **Density by Type**:
  - Iron: 7300-7900 kg/m³
  - Stony-iron: 4400-5000 kg/m³
  - Ordinary chondrites: 3300-3700 kg/m³
  - Carbonaceous chondrites: 2200-3000 kg/m³
- **Typical Entry Velocities**: 11-72 km/s (pre-atmospheric)
- **Terminal Velocities**: 
  - Small meteorites (<10 kg): 50-300 m/s
  - Large meteorites (>1000 kg): 500-3000 m/s
- **Ablation Parameter**: 0.01-0.05 s²/km² (material-dependent)
- **Strength**:
  - Iron: compressive strength 330-430 MPa
  - Stony: compressive strength 150-280 MPa

### 1.2 Impact Energy and Crater Formation

#### Energy Release
- **Kinetic Energy**: E = 0.5mv² [J]
- **Energy-to-TNT Conversion**: 1 kt TNT = 4.184×10¹² J
- **Energy Partitioning**:
  - Shock wave: 50-60%
  - Kinetic energy of ejecta: 30-40%
  - Thermal radiation: 1-3%
  - Seismic waves: 0.5-2%
  - Acoustic waves: <0.1%

#### Airburst vs. Surface Impact
- **Airburst Threshold Condition**: 
  - σ < ρₐv²/2 (where σ is impactor strength, ρₐ is atmospheric density)
- **Altitude of Break-up**: h = H·ln(ρₒρᵢdv²/3σY)
  - Where H is atmospheric scale height (~8 km)
  - ρᵢ is impactor density, d is diameter, σY is yield strength
- **Energy Deposition**: dE/dz = -Cd·ρₐ·A·v³/2
  - Where Cd is drag coefficient, A is cross-sectional area

#### Crater Type by Impactor Characteristics
- **Simple Craters**:
  - Diameter: D < 2-4 km (Earth), D < 15-20 km (Moon)
  - Depth/Diameter ratio: ~1/5
  - Formation time: t ∝ D¹/² (~seconds)
  - Typical for: Small, low-velocity impactors
- **Complex Craters**:
  - Diameter: D > 2-4 km (Earth), D > 15-20 km (Moon)
  - Depth/Diameter ratio: ~1/10 to 1/15
  - Central uplift diameter: ~0.22D
  - Formation time: t ∝ D⁴/³ (~minutes)
  - Typical for: Larger, higher-velocity impactors
- **Multi-ring Basins**:
  - Diameter: D > 100 km
  - Ring spacing ratio: ~√2
  - Formation time: t ∝ D² (~hours)
  - Typical for: Very large, high-velocity impactors

#### Crater Morphology by Impactor Type
- **Iron/Metallic Asteroid Impact**:
  - Greater depth
  - More melt production (10-15× more than cometary impact)
  - More projectile preservation
  - Lower ejecta blanket extent
- **Stony Asteroid Impact**:
  - Moderate depth
  - Moderate melt production
  - Moderate ejecta extent
- **Cometary Impact**:
  - Shallower craters
  - Less melt production (due to ice content)
  - Extensive, thinner ejecta blanket
  - Almost no projectile preservation

### 1.3 Advanced Impact Effects

#### Oblique Impacts
- **Crater Morphology Changes with Impact Angle**:
  - Circular craters even at impact angles >30° from vertical
  - Ellipticity becomes significant at angles <15°
- **Ejecta Distribution**:
  - Asymmetric for θ < 45° from horizontal
  - Downrange ejecta enhancement: ∝ cos(θ)
  - "Forbidden zone" uprange for θ < 20°
- **Energy Coupling Efficiency**: η = sin²θ

#### Porous Target Effects
- **Shock Pressure Attenuation**: Enhanced by factor of 2-5
- **Crater Diameter Reduction**: Dᵨ ≈ D(1-φ)^α
  - Where φ is porosity, α ≈ 1/3
- **Energy Partitioning**: More energy to compaction, less to ejecta

#### Water Impacts
- **Transient Cavity**: Depth ≈ 2-3× impactor diameter
- **Rim Wave Height**:
  - Deep water: hᵣ ≈ 0.014R(r/R)^(-5/4)
  - Where R is transient crater radius, r is distance
- **Tsunami Wave Height Decay**: h ∝ r^(-1) in deep water, h ∝ r^(-1/2) in shallow water
- **Run-up Height**: hᵣᵤₙ = 2.831h₀(tan β)^(1/2)
  - Where h₀ is wave height offshore, β is beach slope

## 2. Karstic and Pseudokarstic Processes

### 2.1 Chemical Dissolution Rates

- **Limestone Dissolution Rate**: 
  - R = k(1-c/cₑq)^n [mol/m²/s]
  - Where c is concentration, cₑq is equilibrium concentration
  - k ≈ 10^(-5) mol/m²/s at pH 7, 25°C
  - n ≈ 1-4 depending on flow regime
- **Dissolution Enhancement by CO₂**:
  - PCO₂ effect: R ∝ (PCO₂)^(1/3)
  - Soil CO₂ typically 10-100× atmospheric levels
- **Temperature Dependence**:
  - Arrhenius equation: k = A·exp(-Eₐ/RT)
  - Eₐ ≈ 8.4 kJ/mol for calcite dissolution
- **Annual Surface Lowering Rates**:
  - Temperate climates: 10-100 mm/kyr
  - Tropical climates: 30-300 mm/kyr
  - Alpine environments: 5-40 mm/kyr

### 2.2 Karst Landform Evolution

- **Doline (Sinkhole) Formation**:
  - Critical roof thickness: h < 0.36r (for collapse)
    - Where r is cavity radius
  - Growth rate: dV/dt ∝ Q·(cₑq-c)
    - Where Q is discharge through system
  - Coalescence law: Spacing ∝ Depth^(1/2)
- **Cave Passage Dimensions**:
  - Width/height ratio as function of recharge mode:
    - Vadose (free surface): 0.2-0.5
    - Phreatic (pressurized): 0.7-1.5
  - Passage cross-section: A ∝ Q^(0.8)·S^(-0.4)
    - Where S is slope
- **Cave System Statistics**:
  - Length-depth ratio: L/D ≈ 20-40
  - Network density: 5-15 km/km²
  - Fractal dimension: D ≈ 1.4-1.8
- **Karst Tower Formation**:
  - Spacing ∝ Bedrock thickness
  - Height ∝ Duration^(1/2)·Rainfall
  - Critical density: towers form when >60% of landscape is dolines

### 2.3 Pseudokarst Processes

- **Piping Erosion Rate**:
  - dV/dt = k·τ^a·(τ-τc)·A
  - Where τ is hydraulic shear stress, τc is critical stress
  - a ≈ 1.5, k depends on soil cohesion
- **Suffosion Rate**:
  - dV/dt ∝ Q·i·(i-ic)
  - Where i is hydraulic gradient, ic is critical gradient
  - ic = (Gs-1)/(1+e)·tanφ
- **Thermokarst Development**:
  - Thaw settlement: S = Hi·θi·Ri
    - Where Hi is ice-rich soil thickness, θi is volumetric ice content
    - Ri is consolidation ratio
  - Lateral erosion rate: dL/dt = k·(T-Tf)^a
    - Where T is temperature, Tf is freezing temperature
- **Lava Tube Formation**:
  - Roof thickness stability: h > 0.14w
    - Where w is tube width
  - Normalized length: L/W ≈ 10-50
    - Where W is width

## 3. Badland Formation and Evolution

### 3.1 Badland Initiation

- **Critical Slope for Rill Formation**:
  - Sc = (d/l)·(τc/τ0)
  - Where d is flow depth, l is hillslope length
  - τc is critical shear stress, τ0 is shear stress
- **Drainage Density Threshold**:
  - Dd > 100 km/km² for badland classification
- **Initiation Factors**:
  - Vegetation cover: critical threshold <15%
  - Rainfall intensity: I > 20-30 mm/h (typical)
  - Soil dispersivity: SAR > 10-13

### 3.2 Morphometric Evolution

- **Drainage Density Increase Rate**:
  - dDd/dt = k·(P/E)·(1-V)·S·(D-Dd)
  - Where P is precipitation, E is potential evapotranspiration
  - V is vegetation cover fraction, S is slope
  - D is maximum potential drainage density
- **Gully Headcut Retreat Rate**:
  - dL/dt = k·Q^a·H^b
  - Where Q is discharge, H is headcut height
  - a ≈ 0.5, b ≈ 1-2
- **Sideslope Evolution**:
  - dS/dt = k1·S^2·(1-S/Smax) - k2·S·V
  - Where S is slope, Smax is maximum stable slope
  - V is vegetation cover

### 3.3 Piping and Tunnel Erosion

- **Threshold Hydraulic Gradient**:
  - ic = (γs/γw)·tanφ'
  - Where γs is soil unit weight, γw is water unit weight
  - φ' is effective friction angle
- **Pipe Growth Rate**:
  - dD/dt = k·(τ-τc)·τ^n
  - Where D is pipe diameter, τ is shear stress
  - n ≈ 1.5-2
- **Roof Collapse Threshold**:
  - Critical pipe diameter: Dc = 2c/γ·tanφ
  - Where c is cohesion, γ is unit weight, φ is friction angle

## 4. Soil Creep and Biogenic Transport Processes

### 4.1 Soil Creep Mechanisms

- **Freeze-Thaw Creep Rate**:
  - dz/dt = k·Nft·d·sin(α)
  - Where Nft is annual freeze-thaw cycles
  - d is displacement per cycle, α is slope angle
- **Temperature-Gradient Creep**:
  - dz/dt = k·(dT/dz)·(dT/dt)·sin(α)
  - Where dT/dz is temperature gradient
  - dT/dt is rate of temperature change
- **Moisture-Change Creep**:
  - dz/dt = k·Nw·ΔW·sin(α)
  - Where Nw is wetting/drying cycles
  - ΔW is moisture content change
- **Combined Diffusion-Like Model**:
  - dz/dt = D·∇²z
  - Where D is diffusivity [m²/yr]:
    - Semiarid: 0.001-0.01
    - Temperate: 0.002-0.03
    - Humid: 0.005-0.05

### 4.2 Biogenic Transport

- **Bioturbation Depth**:
  - h = h0·(1-exp(-t/τ))
  - Where h0 is maximum bioturbation depth
  - τ is characteristic timescale
- **Biogenic Transport Rates by Agent**:
  - Earthworms: 1-5 mm/yr surface accretion
  - Ants: 0.1-1 kg/m²/yr soil transport
  - Termites: 0.05-10 kg/m²/yr soil transport
  - Mammals: 0.02-0.5 kg/m²/yr soil transport
- **Tree Throw Impact**:
  - Soil flux: qs = ρ·F·V·d·cos(α)/2
  - Where F is tree fall frequency [trees/m²/yr]
  - V is root ball volume, d is displacement distance
- **Ecosystem-Scale Transport Coefficient**:
  - K ∝ NPP^(0.3-0.5)
  - Where NPP is net primary productivity

## 5. Hydrocompaction and Subsidence Processes

### 5.1 Hydrocompaction

- **Settlement Due to Wetting**:
  - ΔH = H·Cc·log(σ'v/σ'c)
  - Where H is layer thickness, Cc is compression index
  - σ'v is vertical effective stress, σ'c is preconsolidation stress
- **Collapse Potential**:
  - CP = (e0-ef)/(1+e0)
  - Where e0 is initial void ratio, ef is final void ratio
- **Initiation Threshold**:
  - Moisture content: w > wc
  - Critical saturation degree: Sr > 0.6-0.7
- **Rate of Settlement**:
  - dS/dt = Cv·∇²u
  - Where Cv is coefficient of consolidation, u is pore pressure

### 5.2 Subsidence Due to Fluid Extraction

- **Groundwater Withdrawal Subsidence**:
  - ΔH = α·Δp·H
  - Where α is skeletal compressibility [1/Pa]
  - Δp is pore pressure decline [Pa], H is aquifer thickness [m]
- **Aquifer System Compaction Time Constant**:
  - τ = Ssb²/Kv
  - Where Ss is specific storage, b is aquitard half-thickness
  - Kv is vertical hydraulic conductivity
- **Subsidence Bowl Profile**:
  - s(r) = s0·exp(-r²/2a²)
  - Where s0 is maximum subsidence, a is characteristic radius
  - r is distance from center
- **Maximum Subsidence Rate**:
  - ds/dt = Q·Ss·tan(θ)/π
  - Where Q is pumping rate, θ is cone of depression angle

### 5.3 Ground Rupture and Fissuring

- **Critical Differential Subsidence for Fissuring**:
  - Horizontal strain threshold: εh > 0.01-0.1%
- **Fissure Length Relation**:
  - L ∝ ΔH^(0.5-0.7)
  - Where ΔH is differential subsidence
- **Fissure Width**:
  - w = 2·H·εh·tan(45°-φ/2)
  - Where H is depth to bedrock, φ is friction angle
- **Propagation Rate**:
  - dL/dt ∝ dΔH/dt
  - Where ΔH is differential subsidence

## 6. Periglacial Processes and Landforms

### 6.1 Frost Heave and Sorting

- **Frost Heave Magnitude**:
  - h = α·w·df
  - Where α is frost heave coefficient (1.09-1.4)
  - w is water content, df is freeze depth
- **Frost Front Penetration**:
  - z = λ·√t
  - Where λ is Neumann parameter [m/s^(1/2)]
  - t is freezing duration [s]
- **Sorted Circle Diameter**:
  - D ∝ df^(0.5-0.7)
  - Where df is maximum freeze depth
- **Stone Migration Rate**:
  - dx/dt = k·sin(α)·Nft·d^(-1)
  - Where Nft is freeze-thaw cycles, d is stone diameter
  - α is slope angle

### 6.2 Thermokarst Processes

- **Ice Wedge Degradation Rate**:
  - dV/dt = k·(MAAT-MATc)
  - Where MAAT is mean annual air temperature
  - MATc is critical temperature threshold (-2 to -6°C)
- **Thaw Lake Expansion Rate**:
  - dA/dt = k·(Tw-Tf)·(Rw/H)
  - Where Tw is water temperature, Tf is freezing point
  - Rw is wave energy, H is bank height
- **Critical Thaw Subsidence**:
  - S = Hi·(1-(ρi/ρw))
  - Where Hi is ice thickness, ρi is ice density
  - ρw is water density

### 6.3 Patterned Ground Evolution

- **Pattern Wavelength**:
  - λ = 2π·√(EI/k)
  - Where EI is bending stiffness of active layer
  - k is subgrade reaction modulus
- **Active Layer Detachment Rate**:
  - dA/dt = k·S·(P-ET)·(1-Vf)
  - Where S is slope, P is precipitation
  - ET is evapotranspiration, Vf is vegetation factor
- **Geometric Evolution Rate**:
  - dP/dt = k·(P-Peq)·Nft
  - Where P is pattern order parameter
  - Peq is equilibrium state, Nft is freeze-thaw cycles

## 7. Landscape Evolution Time Constants

### 7.1 System Response Timescales

- **Tectonic Forcing**:
  - Orogenic timescale: τo = h/U
  - Where h is mountain height, U is uplift rate
  - Typically 1-10 Ma
- **Relief Development**:
  - τr = L²/K
  - Where L is system length, K is erodibility
  - Typically 0.01-1 Ma for mountainous regions
- **Isostatic Adjustment**:
  - τi = 2πη/ρmgλ
  - Where η is mantle viscosity, λ is flexural parameter
  - Typically 1-5 ka for post-glacial rebound
- **Knickpoint Migration**:
  - τk = L/V
  - Where L is stream length, V is knickpoint velocity
  - Typically 1-100 ka

### 7.2 Climate Change Response

- **Vegetation Response Time**:
  - τv = 0.1-1 ka (biome shift)
  - τv = 0.01-0.1 ka (community composition)
- **Soil Development Time**:
  - τs = 1-100 ka (complete profile development)
  - τs = 0.1-1 ka (A horizon development)
- **Fluvial System Adjustment**:
  - τf = L/K·A^m·S^n
  - Where L is system length, A is drainage area
  - S is slope, K is erodibility coefficient

### 7.3 Landform Relaxation Times

- **Fault Scarp Degradation**:
  - τd = W²/κ
  - Where W is scarp width, κ is diffusivity
  - Typically 1-10 ka for unconsolidated materials
- **Glacial Feature Adjustment**:
  - Moraine degradation: 1-10 ka
  - Cirque headwall retreat: 10-100 ka
- **Crater Degradation**:
  - τc = R²/κ
  - Where R is crater radius, κ is diffusivity
  - Typically 0.1-1 Ma for 1 km crater on Earth

## 8. Scaling Laws in Landscape Evolution

### 8.1 Spatial Scaling

- **Drainage Network Properties**:
  - Hack's Law: L ∝ A^h, h ≈ 0.5-0.6
  - Horton's Laws:
    - Bifurcation ratio: 3-5
    - Length ratio: 1.5-3.5
    - Area ratio: 3-6
- **Valley Network Properties**:
  - Width-to-depth ratio: w/d ∝ A^m
  - m ≈ 0.3-0.5 for bedrock valleys
  - m ≈ 0.1-0.3 for alluvial valleys
- **Mountain Range Geometry**:
  - Width-to-height ratio: W/H ≈ 5-10
  - Spacing-to-height ratio: S/H ≈ 10-20

### 8.2 Temporal Scaling

- **Erosion Event Frequency-Magnitude**:
  - N(>E) ∝ E^(-β), β ≈ 1-1.5
  - Where N(>E) is number of events exceeding magnitude E
- **Landslide Frequency-Area Distribution**:
  - N(>A) ∝ A^(-α), α ≈ 1.3-1.6
  - Where N(>A) is number of landslides exceeding area A
- **Channel Adjustment Timescale**:
  - τ ∝ L^p
  - Where p ≈ 0.7-1.5 depending on dominant process

### 8.3 Feature Self-Organization

- **Drainage Network Fractal Dimension**:
  - D ≈ 1.8-2.0 (space-filling)
- **Coastline Fractal Dimension**:
  - D ≈ 1.2-1.4
- **Mountain Range Profile Fractal Dimension**:
  - D ≈ 1.2-1.5
- **Topographic Spectral Scaling**:
  - Power spectrum: S(k) ∝ k^(-β)
  - β ≈ 2 (Brownian surface)
  - β ≈ 3 (filtered Brownian surface)