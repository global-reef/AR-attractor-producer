maturity_lookup <- tibble::tribble(
  ~Sci_Name,                     ~Lmax_cm, ~Lmat_cm, ~source, ~notes,
  
  # ----- GROUPERS / SNAPPERS / EMPERORS -----
  "Epinephelus fasciatus",        52,       39.1,     "fishbase",       "Estimated from Linf",
  "Epinephelus fuscoguttatus",    120,      51.8,     "fishbase",       "Estimated from Linf",
  "Lutjanus argentimaculatus",    150,      57,       "fishbase",       "Estimated from Linf",
  "Lutjanus russellii",           50,       29,       "fishbase",       "Estimated from Linf",
  "Lutjanus vitta",               25.5,     15.3,     "fishbase",       "Estimated from Linf",
  "Plectorhinchus chaetodonoides",72,       47.4,     "fishbase",       "Estimated from Linf",
  "Diagramma pictum",             64.3,     31,       "fishbase",       "Measured: off Abu Dhabi, 2000-2003",
  
  # ----- PREDATORS / PELAGICS -----
  "Gnathanodon speciosus",        120,       32.5,    "fishbase",       "Ras El-Kheima, Umm Al-Qwain, Ajman and Sharjah/2014",
  "Sphyraena barracuda",          200,       66,      "fishbase",       "Florida Keys/1996-1998",
  "Sphyraena flavicauda",         60,        28,      "fishbase",       "off Alexandria/1998-1999",
  "Sphyraena jello",              150,       81.7,    "fishbase",       "Estimated from Linf",
  "Sphyraena qenie",              170,       85.6,    "fishbase",       "Estimated from Linf",
  
  # ----- RABBITFISH / EMPERORFISH / PARROTFISH -----
  "Siganus javus",                53,       30.6,     "fishbase",       "Estimated from Linf",
  "Siganus virgatus",             30,       18.5,     "fishbase",       "Estimated from Linf",
  "Scarus rivulatus",             40,       21.8,     "fishbase",       "Estimated from Linf",
  
  # ----- BUTTERFLYFISH / ANGELFISH / DAMSELFISH -----
  "Chaetodon octofasciatus",      12,       8.2,      "fishbase",       "Estimated from Linf",
  "Chaetodon weibeli",            19,       12.4,     "fishbase",       "Estimated from Linf",
  "Heniochus acuminatus",         25,       15.7,     "fishbase",       "Estimated from Linf",
  "Pomacanthus annularis",        45,       26.5,     "fishbase",       "Estimated from Linf",
  "Pomacentrus alexanderae",      9,        6.5,      "fishbase",       "Estimated from Linf",
  "Neopomacentrus cyanomos",      10,       7,        "fishbase",       "Estimated from Linf",
  
  # ----- WRASSES / CLEANERS -----
  "Thalassoma lunare",            45,       26.5,     "fishbase",       "Estimated from Linf",
  "Cheilinus fasciatus",          40,       23.8,     "fishbase",       "Estimated from Linf",
  "Labroides dimidiatus",         14,       9.4,      "fishbase",       "Estimated from Linf",
  
  # ----- TRIGGERFISH / FUSILIERS -----
  "Ballistoides viridescens",     75,       41.5,     "fishbase",       "Estimated from Linf",
  "Caesio xanthonota",            40,       23.8,     "fishbase",       "Estimated from Linf",
  
  # ----- OTHER / MIXED -----
  "Flavocaranx bajad",            55,       24.8,     "fishbase",       "Estimated from Linf", # old: Carangoides bajad
  "Turrum fulvoguttatum",         120,      69.2,     "fishbase",       "Estimated from Linf", # old: Carangoides fulvoguttatum
  "Plectorhinchus gibbosus",      75,       45.9,     "fishbase",       "Estimated from Linf", 
  "Plectropomus spp.",            120,      37.3,     "fishbase",       "Used P. leopardus, from GBR",
  "Ephippidae spp.",              70,       39.1,     "fishbase",       "Used P. teira, Estimated from Linf"
)

