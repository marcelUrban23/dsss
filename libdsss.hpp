#ifndef DSSS_HPP
#define DSSS_HPP

#include <vector>
#include <numeric>
#include <complex>
#include <algorithm>
#include <iostream>

#include "../01_libfft/fft2.hpp"


class Dsss
{
private:
    fft2<std::complex<float>> *m_fft;

    uint64_t  m_fft_leng,
              m_fft_exp;

public:
    Dsss()
    {
        m_fft_exp  = 16;
        m_fft_leng = 1 << m_fft_exp;

        m_fft = new fft2<std::complex<float>>(m_fft_exp);
    }


    /// @brief Berechnet die Autokorrelierte des Eingangsvektors
    std::vector<std::complex<float>>
    getAutoCorr(std::vector<std::complex<float>> input)
    {
        std::vector<std::complex<float>> output(m_fft->fft(input.data()));

        std::transform(output.begin(), output.end(), output.begin(),
                       [](std::complex<float> c)
                       {return std::conj(c) * c;});

        return m_fft->fft(output.data(), 1);
    }



    /// @brief Kreuzkorrelierte zweier Signale berechnen.
    std::vector<std::complex<float>>
    getCrossCorr(std::vector<std::complex<float>> input,
                 std::vector<std::complex<float>> input_2)
    {
        std::vector<std::complex<float>> output(m_fft->fft(input.data()));

        std::transform(output.begin(), output.end(),
                       m_fft->fft(input_2.data()).begin(),
                       output.begin(),
                       [](std::complex<float> val, std::complex<float> val_2)
                       {return (val * std::conj(val_2));});


        return m_fft->fft(output.data(), 1);
    }



    /// @brief Berechnet die Absolutwerte des Eingangsvectors
    std::vector<double>
    getAbs(std::vector<std::complex<float>> input)
    {
        std::vector<double> output(input.size());

        std::transform(input.begin(), input.end(), output.begin(),
                       [](std::complex<float> c)
                       {return std::abs(c);});

        return output;
    }



    /// @brief Berechnet das Leistungsdichtespektrum des Eingangsvectors
    ///        == Sxx = fft(x) * conj(fft(x))
    std::vector<std::complex<float>>
    getPSD(std::vector<std::complex<float>> input)
    {
        std::vector<std::complex<float>> output(m_fft->fft(input.data()));

        std::transform(output.begin(), output.end(), output.begin(),
                       [](std::complex<float> c)
                       {return std::conj(c) * c;});

        return output;
    }



    /// @brief Versuvht die Chiplaenge anhand der Autokorrelierten zu berechnen
    /// @param input Absolutwerte der Autokorreliertens
    /// @return statistisch bester Treffer
    uint64_t
    estimateChipLeng(std::vector<double> input, double threshold = 5.0)
    {
        std::vector<uint64_t> tmp_0,
                              tmp_1;

        tmp_0.reserve(input.size());
        tmp_1.reserve(input.size());

        uint64_t guard_interval = 10,
                 min_chip_leng  = 20,
                 max_chip_leng  = 65536,

                 w              =  0,
                 x              =  0;

        double average = std::accumulate(input.begin(), input.end(), 0.0)
                       / input.size();


        for(w = guard_interval; w < (input.size() - guard_interval); ++w)
        {
            if(   (input[w] > input[w - 1])
               && (input[w] > input[w + 1])
               && (input[w] > (average * threshold)))
            {
                tmp_0.push_back(w - x);

                x = w;
            }
        }

        if(tmp_0.empty())
        {
#ifdef DEBUG
            std::cerr << "estimateChipLeng(): Fehlanzeige" << std::endl;
#endif
            return 0;
        }


        std::sort(tmp_0.begin(),
                  tmp_0.end(),
                  std::less<uint64_t>());

        for(auto u : tmp_0)
        {
            int is_new = 1;
            for(auto v : tmp_1)
            {
                if(u == v)
                {
                    is_new = 0;
                }
            }
            if(   is_new
               && (u > min_chip_leng)
               && (u < max_chip_leng))
            {
                tmp_1.push_back(u);
            }
        }

        std::vector<uint64_t> tmp_2(tmp_1.size());

        for(w = 0; w < tmp_0.size(); ++w)
        {
            for(x = 0; x < tmp_1.size(); ++x)
            {
                if(0 == (tmp_0[w] % tmp_1[x]))
                {
                    tmp_2[x]++;
                    break;
                }
            }
        }

        std::vector<uint64_t>::iterator best = std::max_element(tmp_2.begin(),
                                                                tmp_2.end());

        return static_cast<int>(tmp_1[ static_cast<uint64_t>(
                                         std::distance(tmp_2.begin(), best))]);
    }



    /// @brief Bestimmt die Polaritaet des Chipmusters
    /// @param IQ Eingangssignal
    /// @param Chiplaenge
    /// @return Chipmuster
    std::vector<int>
    estimateChipPattern(std::vector<std::complex<float>> input,
                        uint64_t                         chip_leng)
    {
        uint64_t chip_cnt = input.size() / chip_leng;

        std::vector<double> phases(chip_leng, 0.0);

        for(uint64_t w = 0; w < chip_cnt; ++w)
        {
            float phi = 0.0;
            for(uint64_t x = 0; x < chip_leng; ++x)
            {
                phi += std::arg(input[w * chip_leng + x]);
            }
            phi /= chip_leng;

            for(uint64_t x = 0; x < chip_leng; ++x)
            {
                phases[x] += static_cast<double>((std::arg(input[  w
                                                                 * chip_leng
                                                                 + x]) - phi));
            }
        }

        std::vector<int> output(chip_leng, 0);

        for(uint64_t w = 0; w < chip_leng; ++w)
        {
            (phases[w] / chip_cnt) < 0.0 ? (output[w] = -1)
                                         : (output[w] = 1);
        }


        return output;
    }



    /// @brief Versucht eine dominante Spitze zu finden
    /// @param input Absolutwert der Eingangsdaten
    /// @param guard_interval Abstand zum Anfang und Ende
    /// @param peak_to_side_ratio Adresse fuer Verhaeltnis
    uint64_t
    estimateQualifiedPeak(std::vector<double> input,
                          int64_t             guard_interval,
                          double*             peak_to_side_ratio)
    {
        uint64_t space_to_peak = 5,
                 area_to_peak  = 25;

        std::vector<double>::iterator iter =
                std::max_element(input.begin() + guard_interval,
                                 input.end()   - guard_interval);

        double peak = *iter;

        uint64_t peak_pos = static_cast<uint64_t>(std::distance(input.begin(),
                                                                iter));

        double side_average = 0.0;
        for(uint64_t w = 0; w < area_to_peak; ++w)
        {
            side_average += (  input[peak_pos - space_to_peak - w]
                             + input[peak_pos - space_to_peak + w]);
        }
        side_average /= (2 * area_to_peak);

        peak_to_side_ratio[0] = peak / side_average;



        return peak_pos;
    }


    /// @brief Enspreizt input mittels der chip Sequent
    /// @param input Signal Abtastwerte
    /// @param chip Spreiz- / Entspreizmuster [1,-1]
    /// @return entspreizte Abtastwerte
    ///
    std::vector<std::complex<float>>
    despread(std::vector<std::complex<float>> input, std::vector<int> chip)
    {
        std::vector<std::complex<float>> fc_tmp(input.size(),
                                                std::complex<float>(0.0, 0.0));

        std::transform(chip.begin(), chip.end(), fc_tmp.begin(),
                       [](int val) -> std::complex<float>
                       {return val;});

        double ratio = 0.0;

        uint64_t pos = estimateQualifiedPeak(getAbs(getCrossCorr(input, fc_tmp)),
                                             static_cast<int64_t>(chip.size()),
                                             &ratio);

        uint64_t cnt = (  chip.size()
                          - (static_cast<uint64_t>(pos) % chip.size()))
                       % chip.size();


        std::vector<std::complex<float>> output(input.size());

        for(uint64_t w = 0; w < input.size(); ++w)
        {
            if(chip[cnt] < 0)
            {
                output[w].real(-input[w].real());
                output[w].imag(-input[w].imag());
            }
            else
            {
                output[w] = input[w];
            }

            if(++cnt > (chip.size() - 1))
            {
                cnt = 0;
            }
        }


        return output;
    }
};

#endif




