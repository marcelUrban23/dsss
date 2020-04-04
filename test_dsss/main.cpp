#include <QCoreApplication>

#include <QDataStream>

#include <QFile>

#include <iostream>

#include <stdio.h>

#include "../libdsss.hpp"

void
foo(void);

void
CCtoFC(QString file_in_name);


int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    foo();

    //CCtoFC("/home/hubert/01_PROJECTS/14_DSSS/test_dsss/hackrf.sc");


    return 0;
}




void
foo(void)
{
    uint64_t bytes_to_read = 1 << 16;

    Dsss *dsss = new Dsss();


//    QFile file_in;
//    file_in.open(stdin, QIODevice::ReadOnly);

//    QFile file_out(
//        "/home/hubert/01_PROJECTS/14_DSSS/test_dsss/test_rec_101m88_10m0.32fc");
//    file_out.open(QIODevice::WriteOnly);


    std::vector<char> input(bytes_to_read);
    std::vector<std::complex<float>> tmp(input.size() / 2);

    uint64_t bytes_read = fread(reinterpret_cast<char*>(input.data()),
                                sizeof(char),
                                bytes_to_read,
                                stdin);


    if(bytes_to_read != bytes_read)
    {
        std::cerr << "FEHLER readData(): " << bytes_read << " / " <<
                     bytes_to_read << std::endl;
    }
    else
    {
        std::cerr << "bytes_read: " << bytes_read << std::endl;
    }


    for(uint64_t x = 0; x < tmp.size(); ++x)
    {
        tmp[x].real(static_cast<float>(input[2 * x    ]));
        tmp[x].imag(static_cast<float>(input[2 * x + 1]));
    }

//    file_out.write(reinterpret_cast<const char*>(tmp.data()),
//                   static_cast<int64_t>(  tmp.size()
//                                        * sizeof(std::complex<float>)));


     uint64_t i = dsss->estimateChipLeng(dsss->getAbs(tmp));

     if(i > 0)
     {
         std::cerr << "i: " << i;
         std::vector<int> chip(dsss->estimateChipPattern(tmp, i));

         std::vector<std::complex<float>> vcf_tmp(input.size());

         std::transform(chip.begin(), chip.end(), vcf_tmp.begin(),
                        [](int val) -> std::complex<float>
                        {return val;});

         double ratio = 0.0;

         uint64_t pos = dsss->estimateQualifiedPeak(
                            dsss->getAbs(
                                dsss->getCrossCorr(tmp,
                                                   vcf_tmp)),
                                                   static_cast<int64_t>(
                                                                   chip.size()),
                                                   &ratio);

         std::cerr << "---> ratio: " << ratio << " (" << pos << "}" <<
         std::endl;
         exit(EXIT_SUCCESS);
     }


    return;
}



/// @brief Umwandlung von interleaved char (ichar) zu std::complex<float>
void
CCtoFC(QString file_in_name)
{
    int64_t bytes_to_read = 1 << 16;

    std::vector<char> input(static_cast<uint64_t>(bytes_to_read));


    QFile file_in(file_in_name);
    file_in.open(QIODevice::ReadOnly);

    QFile file_out("out.32fc");
    file_out.open(QIODevice::WriteOnly);

    std::vector<std::complex<float>> tmp(input.size() / 2);

    while(static_cast<int>(bytes_to_read) == file_in.read(input.data(),
                                                          bytes_to_read))
    {
        for(uint64_t x = 0; x < tmp.size(); ++x)
        {
            tmp[x].real(static_cast<float>(input[2 * x    ]));
            tmp[x].imag(static_cast<float>(input[2 * x + 1]));
        }

        file_out.write(reinterpret_cast<const char*>(tmp.data()),
                       static_cast<int64_t>(  tmp.size()
                                            * sizeof(std::complex<float>)));
    }

    file_in.close();
    file_out.close();


    std::cerr << "end" << std::endl;
    exit(EXIT_SUCCESS);
}
