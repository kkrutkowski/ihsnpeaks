#include <QApplication>
#include <QImage>
#include <QPainter>
#include <cstdio>
#include <cmath>
#include <vector>

#include "windows/phased_lightcurve.h"
#include "lc_period.h"

static double calc_median(const float *y, uint32_t n) {
    if (n == 0) return 0.0;
    std::vector<float> ys(y, y + n);
    std::sort(ys.begin(), ys.end());
    return (n % 2 == 1) ? ys[n / 2] : 0.5 * (ys[n / 2 - 1] + ys[n / 2]);
}

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    // Test 1: BLAP-035 (pulsating star)
    {
        const char *filename = "test_data/OGLE-BLAP-035.dat";
        double freq = 49.69881;
        lc_data_t d;
        if (lc_load_dat(filename, &d) == 0 && d.n > 0) {
            lc_detrend(&d);
            lc_periodogram_config_t cfg;
            memset(&cfg, 0, sizeof(cfg));
            cfg.method = LC_SPEC_AOV;
            cfg.nterms = 3;
            cfg.oversampling = 5.0;

            std::vector<float> model(d.n);
            lc_model_style_t style = LC_MODEL_LINE;
            lc_compute_phased_model(&d, &cfg, freq, model.data(), &style);

            // Phase offset logic:
            // d_bright vs d_faint
            double med = calc_median(d.y, d.n);
            uint32_t min_idx = 0, max_idx = 0;
            float min_val = model[0], max_val = model[0];
            for (uint32_t i = 0; i < d.n; ++i) {
                if (model[i] < min_val) { min_val = model[i]; min_idx = i; }
                if (model[i] > max_val) { max_val = model[i]; max_idx = i; }
            }
            double d_bright = std::abs((double)min_val - med) * LC_BRIGHTNESS_BIAS;
            double d_faint  = std::abs((double)max_val - med) * 1.0;
            // Pulsating (d_bright >= d_faint): trough (max_idx) at 0.0
            // Eclipsing (d_faint > d_bright): trough (max_idx) at 0.5
            double target_phase = (d_faint > d_bright) ? 0.5 : 0.0;
            uint32_t target_idx = max_idx;

            double t0 = (d.x[0] + d.x[d.n - 1]) / 2.0;
            double phi_target = std::fmod((d.x[target_idx] - t0) * freq, 1.0);
            if (phi_target < 0.0) phi_target += 1.0;
            double offset = std::fmod(target_phase - phi_target, 1.0);
            if (offset < 0.0) offset += 1.0;

            printf("BLAP-035 (Pulsating) offset = %f, trough target = %f\n", offset, target_phase);

            PhasedLightCurveWidget w("Phased Light Curve - BLAP AoV");
            w.resize(400, 300);
            w.setData(&d);
            w.setFrequency(freq);
            w.setPhaseOffset(offset);
            w.setModel(model.data(), d.n, style, freq);
            w.setDisplayModel(true);

            QImage img(400, 300, QImage::Format_ARGB32_Premultiplied);
            img.fill(QColor(26, 26, 31));
            QPainter p(&img);
            w.render(&p);
            p.end();
            img.save("/home/krutkowski/Pulpit/ihsnpeaks-dev/ihsnpeaks/downstream/lc-qt/test_phased_verify/blap_trough0.png");
            lc_free(&d);
        }
    }

    // Test 2: ECL-04969 (eclipsing binary)
    {
        const char *filename = "test_data/OGLE-LMC-ECL-04969.dat";
        double freq = 0.173363;
        lc_data_t d;
        if (lc_load_dat(filename, &d) == 0 && d.n > 0) {
            lc_detrend(&d);
            lc_periodogram_config_t cfg;
            memset(&cfg, 0, sizeof(cfg));
            cfg.method = LC_SPEC_AOV;
            cfg.nterms = 3;
            cfg.oversampling = 5.0;

            std::vector<float> model(d.n);
            lc_model_style_t style = LC_MODEL_LINE;
            lc_compute_phased_model(&d, &cfg, freq, model.data(), &style);

            double med = calc_median(d.y, d.n);
            uint32_t min_idx = 0, max_idx = 0;
            float min_val = model[0], max_val = model[0];
            for (uint32_t i = 0; i < d.n; ++i) {
                if (model[i] < min_val) { min_val = model[i]; min_idx = i; }
                if (model[i] > max_val) { max_val = model[i]; max_idx = i; }
            }
            double d_bright = std::abs((double)min_val - med) * LC_BRIGHTNESS_BIAS;
            double d_faint  = std::abs((double)max_val - med) * 1.0;
            // Eclipsing: trough (max_idx) at 0.5
            double target_phase = (d_faint > d_bright) ? 0.5 : 0.0;
            uint32_t target_idx = max_idx;

            double t0 = (d.x[0] + d.x[d.n - 1]) / 2.0;
            double phi_target = std::fmod((d.x[target_idx] - t0) * freq, 1.0);
            if (phi_target < 0.0) phi_target += 1.0;
            double offset = std::fmod(target_phase - phi_target, 1.0);
            if (offset < 0.0) offset += 1.0;

            printf("ECL-04969 (Eclipsing) offset = %f, trough target = %f\n", offset, target_phase);

            PhasedLightCurveWidget w("Phased Light Curve - ECL AoV");
            w.resize(400, 300);
            w.setData(&d);
            w.setFrequency(freq);
            w.setPhaseOffset(offset);
            w.setModel(model.data(), d.n, style, freq);
            w.setDisplayModel(true);

            QImage img(400, 300, QImage::Format_ARGB32_Premultiplied);
            img.fill(QColor(26, 26, 31));
            QPainter p(&img);
            w.render(&p);
            p.end();
            img.save("/home/krutkowski/Pulpit/ihsnpeaks-dev/ihsnpeaks/downstream/lc-qt/test_phased_verify/ecl_trough05.png");
            lc_free(&d);
        }
    }

    return 0;
}
