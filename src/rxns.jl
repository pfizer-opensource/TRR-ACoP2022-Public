function get_rxns(df)
       rn = @reaction_network LipoModel begin
        @species begin
            a100(t) = $(df.TV[df.Parameter .== "a100"][1])
            tg_a100(t) = $(df.TV[df.Parameter .== "tg_a100"][1])
            ch_a100(t) = $(df.TV[df.Parameter .== "ch_a100"][1])
            a1(t) = $(df.TV[df.Parameter .== "a1"][1])
            tg_a1(t) = $(df.TV[df.Parameter .== "tg_a1"][1])
            ch_a1(t) = $(df.TV[df.Parameter .== "ch_a1"][1])
            a48(t) = $(df.TV[df.Parameter .== "a48"][1])
            tg_a48(t) = $(df.TV[df.Parameter .== "tg_a48"][1])
            ch_a48(t) = $(df.TV[df.Parameter .== "ch_a48"][1])
            fa_lcy(t) = $(df.TV[df.Parameter .== "fa_lcy"][1])
            fa_ler(t) = $(df.TV[df.Parameter .== "fa_ler"][1])
            tg_lcy(t) = $(df.TV[df.Parameter .== "tg_lcy"][1])
            tg_ler(t) = $(df.TV[df.Parameter .== "tg_ler"][1])
            ch_l(t) = $(df.TV[df.Parameter .== "ch_l"][1])
            ldlr(t) = $(df.TV[df.Parameter .== "ldlr"][1])
            pk9(t) = $(df.TV[df.Parameter .== "pk9"][1])
            ctpi_gut(t) = $(df.TV[df.Parameter .== "ctpi_gut"][1])
            ctpi_cent(t) = $(df.TV[df.Parameter .== "ctpi_cent"][1])
            ctpi_Q1(t) = $( df.TV[df.Parameter .== "ctpi_Q1"][1])
        end
        @parameters begin
            ks_a100 = $(df.TV[df.Parameter .== "ks_a100"][1])
            ks_a1 = $(df.TV[df.Parameter .== "ks_a1"][1])
            kcl_a1 = $(df.TV[df.Parameter .== "kcl_a1"][1])
            krct = $(df.TV[df.Parameter .== "krct"][1])
            ksrb1 = $(df.TV[df.Parameter .== "ksrb1"][1])
            klpl = $(df.TV[df.Parameter .== "klpl"][1])
            kldlr_a100 = $(df.TV[df.Parameter .== "kldlr_a100"][1])
            kldlr_a48 = $(df.TV[df.Parameter .== "kldlr_a48"][1])
            ks_a48 = $(df.TV[df.Parameter .== "ks_a48"][1])
            kctp = $(df.TV[df.Parameter .== "kctp"][1])
            ks_fa = $(df.TV[df.Parameter .== "ks_fa"][1])
            kdnl = $(df.TV[df.Parameter .== "kdnl"][1])
            kd_fa = $(df.TV[df.Parameter .== "kd_fa"][1])
            kest_lcy = $(df.TV[df.Parameter .== "kest_lcy"][1])
            kest_ler = $(df.TV[df.Parameter .== "kest_ler"][1])
            klip = $(df.TV[df.Parameter .== "klip"][1])
            ker = $(df.TV[df.Parameter .== "ker"][1])
            ks_ch = $(df.TV[df.Parameter .== "ks_ch"][1])
            kd_ch = $(df.TV[df.Parameter .== "kd_ch"][1])
            ks_ldlr = $(df.TV[df.Parameter .== "ks_ldlr"][1])
            kd_ldlr = $(df.TV[df.Parameter .== "kd_ldlr"][1])
            ks_pk9 = $(df.TV[df.Parameter .== "ks_pk9"][1])
            kcl_pk9 = $(df.TV[df.Parameter .== "kcl_pk9"][1])
            vd_ler = $(df.TV[df.Parameter .== "vd_ler"][1])
            vd_lcy = $(df.TV[df.Parameter .== "vd_lcy"][1])
            vd_p = $(df.TV[df.Parameter .== "vd_p"][1])
            f_lpl_h = $(df.TV[df.Parameter .== "f_lpl_h"][1])
            alpha_a48 = $(df.TV[df.Parameter .== "alpha_a48"][1])
            beta_a48 = $(df.TV[df.Parameter .== "beta_a48"][1])
            alpha_a100 = $(df.TV[df.Parameter .== "alpha_a100"][1])
            beta_a100 = $(df.TV[df.Parameter .== "beta_a100"][1])
            f_ins = $(df.TV[df.Parameter .== "f_ins"][1])
            K_tgler = $(df.TV[df.Parameter .== "K_tgler"][1])
            K_chl = $(df.TV[df.Parameter .== "K_chl"][1])
            cetp_scale = $(df.TV[df.Parameter .== "cetp_scale"][1])
            pk9_scale = $(df.TV[df.Parameter .== "pk9_scale"][1])
            ka_ctpi = $(df.TV[df.Parameter .== "ka_ctpi"][1])
            kel_ctpi = $(df.TV[df.Parameter .== "kel_ctpi"][1])
            vd_ctpi_cent = $(df.TV[df.Parameter .== "vd_ctpi_cent"][1])
            vd_ctpi_Q1 = $(df.TV[df.Parameter .== "vd_ctpi_Q1"][1])
            Q_ctpi = $(df.TV[df.Parameter .== "Q_ctpi"][1])
        end
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

    end 
    return rn

end

# Design notes:
#   * We do not distinguish between apoB100 particles (e.g., VLDL, IDL, LDL)
#   * We do not distinguish between FC and CE on any particles or in liver
#   * Cholesterol in the liver fuses between a lot of compartments (membranes), we don't attempt to compartmentalize it
#   * Assumed an average of 5 ApoA1s/HDL particle: https://doi.org/10.1016/j.jlr.2021.100099
#   * For readability, a factor of 5 has been absorbed into many of the rate constants that are *[ApoA1]