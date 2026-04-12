/*
 * Deterministic map-number palette for reaction depiction.
 */
#pragma once
#ifndef SMSD_COLOR_PALETTE_HPP
#define SMSD_COLOR_PALETTE_HPP

#include <array>
#include <cstdint>

namespace smsd {

struct PaletteColor {
    std::uint8_t r = 0;
    std::uint8_t g = 0;
    std::uint8_t b = 0;
};

inline PaletteColor mapNumberToPaletteColor(int mapNumber) {
    static constexpr std::array<PaletteColor, 20> palette = {{
        {173, 216, 230}, // light blue
        {122, 201, 111}, // green
        {222, 166,  86}, // amber
        {150, 150, 150}, // gray
        {255, 202,  58}, // gold
        { 79, 172, 232}, // blue
        {214,  93,  93}, // red
        { 62, 100, 180}, // cobalt
        { 93, 174, 222}, // sky
        { 58, 119, 175}, // denim
        { 56,  95, 176}, // royal blue
        {243,  88,  73}, // vermilion
        {164, 107, 199}, // purple
        {151, 174,  84}, // olive
        {205,  92,  92}, // indian red
        {232, 111,  91}, // coral
        {209, 209, 209}, // light gray
        {255,  53,  53}, // strong red
        {255, 140, 105}, // salmon
        { 75, 135, 185}  // steel blue
    }};
    if (mapNumber <= 0) {
        return {21, 101, 192};
    }
    return palette[static_cast<size_t>((mapNumber - 1) % palette.size())];
}

} // namespace smsd

#endif // SMSD_COLOR_PALETTE_HPP
