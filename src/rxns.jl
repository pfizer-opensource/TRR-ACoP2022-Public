function get_rxns()
    rn = @reaction_network LipoModel begin

        # Plasma Compartment Reactions:

        # ApoB100 Particle:
        ks_a100, ∅ --> a100                                     # [mM/day], liver prod. apoB100
        ks_a100 * alpha_a100 * tg_ler / (K_tgler + tg_ler), tg_ler ⇒ vd_ler / vd_p * tg_a100   # [mM/day], liver prod. TG  (apoB100)
        ks_a100 * beta_a100 * ch_l / (K_chl + ch_l), ch_l ⇒ (vd_ler + vd_lcy) / vd_p * ch_a100     # [mM/day], liver prod. CE (apoB100)
        ksrb1, ch_a100 --> vd_p/(vd_ler + vd_lcy) * ch_l                               # [mM/day], SRB1 activity apoB100
        klpl*(1-f_lpl_h), tg_a100 --> ∅                                # [mM/day], lipase activity apoB100
        kldlr_a100 * ldlr, a100 --> ∅                                  # [mM/day], LDL-R activity apoB100
        kldlr_a100 * ldlr, tg_a100 --> 3 * vd_p / vd_lcy * fa_lcy            # [mM/day], LDL-R activity apoB100
        kldlr_a100 * ldlr, ch_a100 --> vd_p / (vd_ler + vd_lcy) * ch_l       # [mM/day], LDL-R activity apoB100

        # ApoA1 particle:
        (ks_a1, kcl_a1), ∅ <--> a1                           # [mM/day], birth/death apoA1
        krct * (a1), ∅ --> ch_a1                                  # [mM/day], RCT CE onto apoA1
        ksrb1, ch_a1 --> vd_p/(vd_ler + vd_lcy) * ch_l                                  # [mM/day], SRB1 activity apoA1
        klpl*(1-f_lpl_h), tg_a1 --> ∅                                  # [mM/day], lipase activity apoB48

        # ApoB48 Particle:
        (ks_a48, ks_a48 * alpha_a48, ks_a48 * beta_a48), ∅ --> (a48, tg_a48, ch_a48) # Birth of apoB48 particles (diet)
        ksrb1, ch_a48 --> vd_p/(vd_ler + vd_lcy) * ch_l                 # [mM/day], SRB1 activity apoB48
        klpl*(1-f_lpl_h)*f_ins, tg_a48 --> ∅                            # [mM/day], lipase activity apo48
        kldlr_a48 * ldlr, a48 --> ∅                                         # [mM/day], LDL-R activity apoB48
        kldlr_a48 * ldlr, tg_a48 --> 3*vd_p / vd_lcy * fa_lcy               # [mM/day], LDL-R activity apoB48
        kldlr_a48 * ldlr, ch_a48 --> vd_p / (vd_ler + vd_lcy) * ch_l        # [mM/day], LDL-R activity apoB48

        # CETP Activity:
        #   *TODO: This section needs a re-think.
        #   *We need to be thinking about ApoB-TG --> ApoA1-TG followed by
        #   *Equimolar ApoA1-Ch --> ApoB-Ch
        kctp*cetp_scale*ch_a1, tg_a100 --> tg_a1      # [mM/day], CETP exchange of TG a100 --> a1
        #kctp*cetp_scale*ch_a1, tg_a48 --> tg_a1        # [mM/day], CETP exchange of TG a48 --> a1
        kctp*cetp_scale*tg_a100, ch_a1 --> ch_a100      # [mM/day], CETP exchange of Ch a1 --> a100
        #kctp*cetp_scale*tg_a1, ch_a1 --> ch_a48        # [mM/day], CETP exchange of Ch a1 --> a48

        # Liver Compartment Reactions:

        # Liver TG:
        ks_fa, ∅ --> fa_lcy # [mM/day], uptake of fatty acids into cytosol
        kdnl, ∅ --> fa_lcy # [mM/day], DNL synthesis of fatty acids
        kd_fa, fa_lcy --> ∅ # [mM/day], beta-oxidation of fatty acids
        kest_lcy, 3 * fa_lcy --> tg_lcy # [mM/day], esterification
        klip, tg_lcy --> 3 * fa_lcy # [mM/day], lipolysis
        ker, fa_lcy --> vd_lcy / vd_ler * fa_ler # [mM/day], uptake into vd_er
        kest_ler, 3*fa_ler --> tg_ler # [mM/day], esterification in ER
        klpl*f_lpl_h, tg_a100 --> 3 * vd_p/vd_lcy*fa_lcy # [mM/day], HL clearance of plasma TG
        klpl*f_lpl_h, tg_a48 --> 3 * vd_p/vd_lcy*fa_lcy # [mM/day], HL clearance of plasma TG
        klpl*f_lpl_h, tg_a1 --> 3 * vd_p/vd_lcy*fa_lcy # [mM/day], HL clearance of plasma TG

        # Liver Cholesterol:
        (ks_ch, kd_ch), ∅ <--> ch_l # [mM/day], birth/death liver cholesterol (HMG-CoA, bile acids)
        (ks_ldlr / ch_l, kd_ldlr * pk9), ∅ <--> ldlr # [mM/day], birth/death LDL-R (SREBP-2, PCSK9-mediated)
        (ks_pk9, kcl_pk9*pk9_scale), ∅ <--> pk9 # [mM/day], birth/death plasma PCSK9

        # Pharmacokinetics:
        ka_ctpi, ctpi_gut --> 1/vd_ctpi_cent * ctpi_cent # [nM/hour], gut absorption 
        kel_ctpi, ctpi_cent --> ∅ # [nM/hour], central clearance
        Q_ctpi/vd_ctpi_cent, ctpi_cent --> vd_ctpi_cent/vd_ctpi_Q1 * ctpi_Q1 # [nM/hour], Central --> Comp1
        Q_ctpi/vd_ctpi_Q1, ctpi_Q1 --> vd_ctpi_Q1/vd_ctpi_cent * ctpi_cent # [nM/hour], Comp1 --> Central

    end ks_a100 ks_a1 kcl_a1 krct ksrb1 klpl kldlr_a100 kldlr_a48 ks_a48 kctp ks_fa kdnl kd_fa kest_lcy kest_ler klip ker ks_ch kd_ch ks_ldlr kd_ldlr ks_pk9 kcl_pk9 vd_ler vd_lcy vd_p alpha_a48 beta_a48 alpha_a100 beta_a100 K_tgler K_chl f_lpl_h f_ins cetp_scale pk9_scale ka_ctpi kel_ctpi vd_ctpi_cent vd_ctpi_Q1 Q_ctpi

    return rn

end

# Design notes:
#   * We do not distinguish between apoB100 particles (e.g., VLDL, IDL, LDL)
#   * We do not distinguish between FC and CE on any particles or in liver
#   * Cholesterol in the liver fuses between a lot of compartments (membranes), we don't attempt to compartmentalize it
#   * Assumed an average of 5 ApoA1s/HDL particle: https://doi.org/10.1016/j.jlr.2021.100099
#   * For readability, a factor of 5 has been absorbed into many of the rate constants that are *[ApoA1]