# Analogy

-   A cathode ray tube creates images by rapidly moving an electron beam
    across a phosphor-coated screen.

-   The beam moves so quickly that it creates the illusion of a stable
    image, even though only one point is illuminated at any given
    instant.

-   If sampled slower than refresh rate, the screen appears as a stable
    image. It may be thought of as the beam being everywhere on the
    screen at once -- a raster superimposition.

-   Increasing sampling rate to near-refresh rate destabilizes the
    apparent image, causing it to flicker, break into bands, and
    ultimately disappear.

-   Matching the sampling window to frame rate, and sampling interval to
    pixel rate, will cause the image to collapse to a single dot --
    raster function collapse.

-   When synchronizing to the dot rate, it is impossible to predict
    where on the screen the dot will appear -- raster uncertainty. But
    once synchronized, sampling with the same timing will always put the
    dot in the same place.