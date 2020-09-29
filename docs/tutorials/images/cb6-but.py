import base64

import cPickle

try:
    from SimpleSession.versions.v65 import (
        beginRestore,
        checkVersion,
        registerAfterModelsCB,
        reportRestoreError,
    )
except ImportError:
    from chimera import UserError

    raise UserError(
        "Cannot open session that was saved in a"
        " newer version of Chimera; update your version"
    )
checkVersion([1, 12, 41613])
from chimera import replyobj

replyobj.status("Restoring session...", blankAfter=0)
replyobj.status("Beginning session restore...", blankAfter=0, secondary=True)
beginRestore()


def restoreCoreModels():
    from SimpleSession.versions.v65 import (
        init,
        restoreColors,
        restoreModelAssociations,
        restoreMolecules,
        restorePseudoBondGroups,
        restoreSurfaces,
        restoreViewer,
        restoreVRML,
    )

    molInfo = cPickle.loads(
        base64.b64decode(
            "gAJ9cQEoVRFyaWJib25JbnNpZGVDb2xvcnECSwFOfYdVCWJhbGxTY2FsZXEDSwFHP9AAAAAAAAB9h1UJcG9pbnRTaXplcQRLAUc/8AAAAAAAAH2HVQVjb2xvcnEFSwFLAH2HVQpyaWJib25UeXBlcQZLAUsAfYdVCnN0aWNrU2NhbGVxB0sBRz/wAAAAAAAAfYdVDG1tQ0lGSGVhZGVyc3EIXXEJTmFVDGFyb21hdGljTW9kZXEKSwFLAX2HVQp2ZHdEZW5zaXR5cQtLAUdAFAAAAAAAAH2HVQZoaWRkZW5xDEsBiX2HVQ1hcm9tYXRpY0NvbG9ycQ1LAU59h1UPcmliYm9uU21vb3RoaW5ncQ5LAUsAfYdVCWF1dG9jaGFpbnEPSwGIfYdVCnBkYlZlcnNpb25xEEsBSwJ9h1UIb3B0aW9uYWxxEX1xElUIb3BlbmVkQXNxE4iJSwEoVQtjYjYtYnV0LnBkYnEUTk5LAXRxFX2Hh3NVD2xvd2VyQ2FzZUNoYWluc3EWSwGJfYdVCWxpbmVXaWR0aHEXSwFHP/AAAAAAAAB9h1UPcmVzaWR1ZUxhYmVsUG9zcRhLAUsAfYdVBG5hbWVxGUsBWAsAAABjYjYtYnV0LnBkYn2HVQ9hcm9tYXRpY0Rpc3BsYXlxGksBiX2HVQ9yaWJib25TdGlmZm5lc3NxG0sBRz/pmZmZmZmafYdVCnBkYkhlYWRlcnNxHF1xHX1xHmFVA2lkc3EfSwFLAEsAhn2HVQ5zdXJmYWNlT3BhY2l0eXEgSwFHv/AAAAAAAAB9h1UQYXJvbWF0aWNMaW5lVHlwZXEhSwFLAn2HVRRyaWJib25IaWRlc01haW5jaGFpbnEiSwGIfYdVB2Rpc3BsYXlxI0sBiH2HdS4="
        )
    )
    resInfo = cPickle.loads(
        base64.b64decode(
            "gAJ9cQEoVQZpbnNlcnRxAksCVQEgfYdVC2ZpbGxEaXNwbGF5cQNLAol9h1UEbmFtZXEESwJYAwAAAENCNn1xBVgDAAAAQlVUXXEGSwFhc4dVBWNoYWlucQdLAlgBAAAAIH2HVQ5yaWJib25EcmF3TW9kZXEISwJLAn2HVQJzc3EJSwKJiYZ9h1UIbW9sZWN1bGVxCksCSwB9h1ULcmliYm9uQ29sb3JxC0sCTn1xDEsBTl1xDUsBSwGGcQ5hhnOHVQVsYWJlbHEPSwJYAAAAAH2HVQpsYWJlbENvbG9ycRBLAk59cRFLAU5dcRJLAUsBhnETYYZzh1UIZmlsbE1vZGVxFEsCSwF9h1UFaXNIZXRxFUsCiX2HVQtsYWJlbE9mZnNldHEWSwJOfYdVCHBvc2l0aW9ucRddcRhLAUsChnEZYVUNcmliYm9uRGlzcGxheXEaSwKJfYdVCG9wdGlvbmFscRt9VQRzc0lkcRxLAkr/////fYd1Lg=="
        )
    )
    atomInfo = cPickle.loads(
        base64.b64decode(
            "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"
        )
    )
    bondInfo = cPickle.loads(
        base64.b64decode(
            "gAJ9cQEoVQVjb2xvcnECS4tOfYdVBWF0b21zcQNdcQQoXXEFKEsDSwRlXXEGKEsDSwVlXXEHKEsDSwZlXXEIKEsDSwdlXXEJKEsESwhlXXEKKEsESwllXXELKEsESwplXXEMKEsFSw5lXXENKEsFSx1lXXEOKEsGSw1lXXEPKEsGSx5lXXEQKEsJSwtlXXERKEsJSw5lXXESKEsKSwxlXXETKEsKSw1lXXEUKEsLSxFlXXEVKEsLSxJlXXEWKEsLS19lXXEXKEsMSxNlXXEYKEsMSxRlXXEZKEsMS2BlXXEaKEsNSw9lXXEbKEsOSxBlXXEcKEsVSxZlXXEdKEsVSxdlXXEeKEsVSxhlXXEfKEsVSxllXXEgKEsWSxplXXEhKEsWSxtlXXEiKEsWSxxlXXEjKEsXSyBlXXEkKEsXSy9lXXElKEsYSx9lXXEmKEsYSzBlXXEnKEsbSx1lXXEoKEsbSyBlXXEpKEscSx5lXXEqKEscSx9lXXErKEsdSyNlXXEsKEsdSyRlXXEtKEseSyVlXXEuKEseSyZlXXEvKEsfSyFlXXEwKEsgSyJlXXExKEsnSyhlXXEyKEsnSyllXXEzKEsnSyplXXE0KEsnSytlXXE1KEsoSyxlXXE2KEsoSy1lXXE3KEsoSy5lXXE4KEspSzJlXXE5KEspS0FlXXE6KEsqSzFlXXE7KEsqS0JlXXE8KEstSy9lXXE9KEstSzJlXXE+KEsuSzBlXXE/KEsuSzFlXXFAKEsvSzVlXXFBKEsvSzZlXXFCKEswSzdlXXFDKEswSzhlXXFEKEsxSzNlXXFFKEsySzRlXXFGKEs5SzplXXFHKEs5SztlXXFIKEs5SzxlXXFJKEs5Sz1lXXFKKEs6Sz5lXXFLKEs6Sz9lXXFMKEs6S0BlXXFNKEs7S0RlXXFOKEs7S1NlXXFPKEs8S0NlXXFQKEs8S1RlXXFRKEs/S0FlXXFSKEs/S0RlXXFTKEtAS0JlXXFUKEtAS0NlXXFVKEtBS0dlXXFWKEtBS0hlXXFXKEtCS0llXXFYKEtCS0plXXFZKEtDS0VlXXFaKEtES0ZlXXFbKEtLS0xlXXFcKEtLS01lXXFdKEtLS05lXXFeKEtLS09lXXFfKEtMS1BlXXFgKEtMS1FlXXFhKEtMS1JlXXFiKEtNS1ZlXXFjKEtNS2VlXXFkKEtOS1VlXXFlKEtOS2ZlXXFmKEtRS1NlXXFnKEtRS1ZlXXFoKEtSS1RlXXFpKEtSS1VlXXFqKEtTS1llXXFrKEtTS1plXXFsKEtUS1tlXXFtKEtUS1xlXXFuKEtVS1dlXXFvKEtWS1hlXXFwKEtdS15lXXFxKEtdS19lXXFyKEtdS2BlXXFzKEtdS2FlXXF0KEteS2JlXXF1KEteS2NlXXF2KEteS2RlXXF3KEtfS2hlXXF4KEtgS2dlXXF5KEtjS2VlXXF6KEtjS2hlXXF7KEtkS2ZlXXF8KEtkS2dlXXF9KEtlS2tlXXF+KEtlS2xlXXF/KEtmS21lXXGAKEtmS25lXXGBKEtnS2llXXGCKEtoS2plXXGDKEtvS3BlXXGEKEtvS3FlXXGFKEtvS3JlXXGGKEtvS3NlXXGHKEtzS3RlXXGIKEtzS3VlXXGJKEtzS3ZlXXGKKEt1S3dlXXGLKEt1S3hlXXGMKEt1S3llXXGNKEt5S3plXXGOKEt5S3tlXXGPKEt5S3xlZVUFbGFiZWxxkEuLWAAAAAB9h1UIaGFsZmJvbmRxkUuLiH2HVQZyYWRpdXNxkkuLRz/JmZmgAAAAfYdVC2xhYmVsT2Zmc2V0cZNLi059h1UIZHJhd01vZGVxlEuLSwF9h1UIb3B0aW9uYWxxlX1VB2Rpc3BsYXlxlkuLSwJ9h3Uu"
        )
    )
    crdInfo = cPickle.loads(
        base64.b64decode(
            "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"
        )
    )
    surfInfo = {
        "category": (1, u"main", {}),
        "probeRadius": (1, 1.4, {}),
        "pointSize": (1, 1, {}),
        "name": [u"MSMS main surface of cb6-but.pdb"],
        "density": (1, 20, {}),
        "colorMode": (1, 1, {}),
        "useLighting": (1, True, {}),
        "transparencyBlendMode": (1, 1, {}),
        "molecule": [0],
        "smoothLines": (1, False, {}),
        "lineWidth": (1, 1, {}),
        "allComponents": (1, True, {}),
        "twoSidedLighting": (1, True, {}),
        "customVisibility": [None],
        "drawMode": (1, 0, {}),
        "display": (1, True, {}),
        "customColors": [(0, None, {})],
    }
    vrmlInfo = {
        "subid": (0, None, {}),
        "display": (0, None, {}),
        "id": (0, None, {}),
        "vrmlString": [],
        "name": (0, None, {}),
    }
    colors = {
        u"Ru": ((0.141176, 0.560784, 0.560784), 1, u"default"),
        u"Re": ((0.14902, 0.490196, 0.670588), 1, u"default"),
        u"Rf": ((0.8, 0, 0.34902), 1, u"default"),
        u"Ra": ((0, 0.490196, 0), 1, u"default"),
        u"Rb": ((0.439216, 0.180392, 0.690196), 1, u"default"),
        u"Rn": ((0.258824, 0.509804, 0.588235), 1, u"default"),
        u"Rh": ((0.0392157, 0.490196, 0.54902), 1, u"default"),
        u"Be": ((0.760784, 1, 0), 1, u"default"),
        u"Ba": ((0, 0.788235, 0), 1, u"default"),
        u"Bh": ((0.878431, 0, 0.219608), 1, u"default"),
        u"Bi": ((0.619608, 0.309804, 0.709804), 1, u"default"),
        u"Bk": ((0.541176, 0.309804, 0.890196), 1, u"default"),
        u"Br": ((0.65098, 0.160784, 0.160784), 1, u"default"),
        u"H": ((1, 1, 1), 1, u"default"),
        u"P": ((1, 0.501961, 0), 1, u"default"),
        u"Os": ((0.14902, 0.4, 0.588235), 1, u"default"),
        u"Es": ((0.701961, 0.121569, 0.831373), 1, u"default"),
        u"Hg": ((0.721569, 0.721569, 0.815686), 1, u"default"),
        u"Ge": ((0.4, 0.560784, 0.560784), 1, u"default"),
        u"Gd": ((0.270588, 1, 0.780392), 1, u"default"),
        u"Ga": ((0.760784, 0.560784, 0.560784), 1, u"default"),
        u"Pr": ((0.85098, 1, 0.780392), 1, u"default"),
        u"Pt": ((0.815686, 0.815686, 0.878431), 1, u"default"),
        u"Pu": ((0, 0.419608, 1), 1, u"default"),
        u"C": ((0.564706, 0.564706, 0.564706), 1, u"default"),
        u"Pb": ((0.341176, 0.34902, 0.380392), 1, u"default"),
        u"Pa": ((0, 0.631373, 1), 1, u"default"),
        u"Pd": ((0, 0.411765, 0.521569), 1, u"default"),
        u"Cd": ((1, 0.85098, 0.560784), 1, u"default"),
        u"Po": ((0.670588, 0.360784, 0), 1, u"default"),
        u"Pm": ((0.639216, 1, 0.780392), 1, u"default"),
        u"Hs": ((0.901961, 0, 0.180392), 1, u"default"),
        u"Ho": ((0, 1, 0.611765), 1, u"default"),
        u"Hf": ((0.301961, 0.760784, 1), 1, u"default"),
        u"K": ((0.560784, 0.25098, 0.831373), 1, u"default"),
        u"He": ((0.85098, 1, 1), 1, u"default"),
        u"Md": ((0.701961, 0.0509804, 0.65098), 1, u"default"),
        u"Mg": ((0.541176, 1, 0), 1, u"default"),
        u"Mo": ((0.329412, 0.709804, 0.709804), 1, u"default"),
        u"Mn": ((0.611765, 0.478431, 0.780392), 1, u"default"),
        u"O": ((1, 0.0509804, 0.0509804), 1, u"default"),
        u"Mt": ((0.921569, 0, 0.14902), 1, u"default"),
        u"S": ((1, 1, 0.188235), 1, u"default"),
        u"W": ((0.129412, 0.580392, 0.839216), 1, u"default"),
        u"Zn": ((0.490196, 0.501961, 0.690196), 1, u"default"),
        u"Eu": ((0.380392, 1, 0.780392), 1, u"default"),
        u"Zr": ((0.580392, 0.878431, 0.878431), 1, u"default"),
        u"Er": ((0, 0.901961, 0.458824), 1, u"default"),
        u"Ni": ((0.313725, 0.815686, 0.313725), 1, u"default"),
        u"No": ((0.741176, 0.0509804, 0.529412), 1, u"default"),
        u"Na": ((0.670588, 0.360784, 0.94902), 1, u"default"),
        u"Nb": ((0.45098, 0.760784, 0.788235), 1, u"default"),
        u"Nd": ((0.780392, 1, 0.780392), 1, u"default"),
        u"Ne": ((0.701961, 0.890196, 0.960784), 1, u"default"),
        u"Np": ((0, 0.501961, 1), 1, u"default"),
        u"Fr": ((0.258824, 0, 0.4), 1, u"default"),
        u"Fe": ((0.878431, 0.4, 0.2), 1, u"default"),
        u"Fm": ((0.701961, 0.121569, 0.729412), 1, u"default"),
        u"B": ((1, 0.709804, 0.709804), 1, u"default"),
        u"F": ((0.564706, 0.878431, 0.313725), 1, u"default"),
        u"Sr": ((0, 1, 0), 1, u"default"),
        u"cornflower blue": ((0.392157, 0.584314, 0.929412), 1, u"default"),
        u"N": ((0.188235, 0.313725, 0.972549), 1, u"default"),
        u"Kr": ((0.360784, 0.721569, 0.819608), 1, u"default"),
        u"Si": ((0.941176, 0.784314, 0.627451), 1, u"default"),
        u"Sn": ((0.4, 0.501961, 0.501961), 1, u"default"),
        u"Sm": ((0.560784, 1, 0.780392), 1, u"default"),
        u"V": ((0.65098, 0.65098, 0.670588), 1, u"default"),
        u"Sc": ((0.901961, 0.901961, 0.901961), 1, u"default"),
        u"Sb": ((0.619608, 0.388235, 0.709804), 1, u"default"),
        u"Sg": ((0.85098, 0, 0.270588), 1, u"default"),
        u"Se": ((1, 0.631373, 0), 1, u"default"),
        u"Co": ((0.941176, 0.564706, 0.627451), 1, u"default"),
        u"Cm": ((0.470588, 0.360784, 0.890196), 1, u"default"),
        u"Cl": ((0.121569, 0.941176, 0.121569), 1, u"default"),
        u"Ca": ((0.239216, 1, 0), 1, u"default"),
        u"Cf": ((0.631373, 0.211765, 0.831373), 1, u"default"),
        u"Ce": ((1, 1, 0.780392), 1, u"default"),
        u"Xe": ((0.258824, 0.619608, 0.690196), 1, u"default"),
        u"Lu": ((0, 0.670588, 0.141176), 1, u"default"),
        u"Cs": ((0.341176, 0.0901961, 0.560784), 1, u"default"),
        u"Cr": ((0.541176, 0.6, 0.780392), 1, u"default"),
        u"Cu": ((0.784314, 0.501961, 0.2), 1, u"default"),
        u"La": ((0.439216, 0.831373, 1), 1, u"default"),
        u"Li": ((0.8, 0.501961, 1), 1, u"default"),
        u"Tl": ((0.65098, 0.329412, 0.301961), 1, u"default"),
        u"Tm": ((0, 0.831373, 0.321569), 1, u"default"),
        u"Lr": ((0.780392, 0, 0.4), 1, u"default"),
        u"Th": ((0, 0.729412, 1), 1, u"default"),
        u"Ti": ((0.74902, 0.760784, 0.780392), 1, u"default"),
        u"Te": ((0.831373, 0.478431, 0), 1, u"default"),
        u"Tb": ((0.188235, 1, 0.780392), 1, u"default"),
        u"Tc": ((0.231373, 0.619608, 0.619608), 1, u"default"),
        u"Ta": ((0.301961, 0.65098, 1), 1, u"default"),
        u"Yb": ((0, 0.74902, 0.219608), 1, u"default"),
        u"Db": ((0.819608, 0, 0.309804), 1, u"default"),
        u"Dy": ((0.121569, 1, 0.780392), 1, u"default"),
        u"I": ((0.580392, 0, 0.580392), 1, u"default"),
        u"U": ((0, 0.560784, 1), 1, u"default"),
        u"Y": ((0.580392, 1, 1), 1, u"default"),
        u"Ac": ((0.439216, 0.670588, 0.980392), 1, u"default"),
        u"Ag": ((0.752941, 0.752941, 0.752941), 1, u"default"),
        u"Ir": ((0.0901961, 0.329412, 0.529412), 1, u"default"),
        u"Am": ((0.329412, 0.360784, 0.94902), 1, u"default"),
        u"Al": ((0.74902, 0.65098, 0.65098), 1, u"default"),
        u"As": ((0.741176, 0.501961, 0.890196), 1, u"default"),
        u"Ar": ((0.501961, 0.819608, 0.890196), 1, u"default"),
        u"Au": ((1, 0.819608, 0.137255), 1, u"default"),
        u"At": ((0.458824, 0.309804, 0.270588), 1, u"default"),
        u"In": ((0.65098, 0.458824, 0.45098), 1, u"default"),
    }
    materials = {u"default": ((0, 0, 0), 30)}
    pbInfo = {
        "category": [u"distance monitor"],
        "bondInfo": [
            {
                "color": (0, None, {}),
                "atoms": [],
                "label": (0, None, {}),
                "halfbond": (0, None, {}),
                "labelColor": (0, None, {}),
                "labelOffset": (0, None, {}),
                "drawMode": (0, None, {}),
                "display": (0, None, {}),
            }
        ],
        "lineType": (1, 2, {}),
        "color": (1, 6, {}),
        "optional": {"fixedLabels": (True, False, (1, False, {}))},
        "display": (1, True, {}),
        "showStubBonds": (1, False, {}),
        "lineWidth": (1, 1, {}),
        "stickScale": (1, 1, {}),
        "id": [-2],
    }
    modelAssociations = {}
    colorInfo = (
        9,
        (u"H", (1, 1, 1, 1)),
        {
            (u"green", (0, 1, 0, 1)): [1],
            (u"", (0.106829, 0.702586, 0.652042, 1)): [0],
            (u"N", (0.188235, 0.313725, 0.972549, 1)): [2],
            (u"", (1, 1, 1, 1)): [7],
            (u"O", (1, 0.0509804, 0.0509804, 1)): [4],
            (u"yellow", (1, 1, 0, 1)): [6],
            (u"cornflower blue", (0.392157, 0.584314, 0.929412, 1)): [5],
            (u"", (0.545455, 0, 1, 1)): [8],
        },
    )
    viewerInfo = {
        "cameraAttrs": {
            "center": (0.031, -0.0945, 1.835),
            "fieldOfView": 17.183565421784,
            "nearFar": (9.0768353811706, -5.4068355123077),
            "ortho": True,
            "eyeSeparation": 50.8,
            "focal": 1.835,
        },
        "viewerAttrs": {
            "silhouetteColor": None,
            "clipping": False,
            "showSilhouette": False,
            "showShadows": False,
            "viewSize": 11.391911504425,
            "labelsOnTop": True,
            "depthCueRange": (0.5, 1),
            "silhouetteWidth": 2,
            "singleLayerTransparency": True,
            "shadowTextureSize": 2048,
            "backgroundImage": [None, 1, 2, 1, 0, 0],
            "backgroundGradient": [
                ("Chimera default", [(1, 1, 1, 1), (0, 0, 1, 1)], 1),
                1,
                0,
                0,
            ],
            "depthCue": True,
            "highlight": 0,
            "scaleFactor": 1.1936205300669,
            "angleDependentTransparency": True,
            "backgroundMethod": 0,
        },
        "viewerHL": 8,
        "cameraMode": "mono",
        "detail": 1.5,
        "viewerFog": None,
        "viewerBG": 7,
    }

    replyobj.status("Initializing session restore...", blankAfter=0, secondary=True)
    from SimpleSession.versions.v65 import expandSummary

    init(dict(enumerate(expandSummary(colorInfo))))
    replyobj.status("Restoring colors...", blankAfter=0, secondary=True)
    restoreColors(colors, materials)
    replyobj.status("Restoring molecules...", blankAfter=0, secondary=True)
    restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo)
    replyobj.status("Restoring surfaces...", blankAfter=0, secondary=True)
    restoreSurfaces(surfInfo)
    replyobj.status("Restoring VRML models...", blankAfter=0, secondary=True)
    restoreVRML(vrmlInfo)
    replyobj.status("Restoring pseudobond groups...", blankAfter=0, secondary=True)
    restorePseudoBondGroups(pbInfo)
    replyobj.status("Restoring model associations...", blankAfter=0, secondary=True)
    restoreModelAssociations(modelAssociations)
    replyobj.status("Restoring camera...", blankAfter=0, secondary=True)
    restoreViewer(viewerInfo)


try:
    restoreCoreModels()
except:
    reportRestoreError("Error restoring core models")

    replyobj.status("Restoring extension info...", blankAfter=0, secondary=True)


try:
    from StructMeasure.DistMonitor import restoreDistances

    registerAfterModelsCB(restoreDistances, 1)
except:
    reportRestoreError("Error restoring distances in session")


def restoreMidasBase():
    formattedPositions = {}
    import Midas

    Midas.restoreMidasBase(formattedPositions)


try:
    restoreMidasBase()
except:
    reportRestoreError("Error restoring Midas base state")


def restoreMidasText():
    from Midas import midas_text

    midas_text.aliases = {}
    midas_text.userSurfCategories = {}


try:
    restoreMidasText()
except:
    reportRestoreError("Error restoring Midas text state")


def restore_cap_attributes():
    cap_attributes = {
        "cap_attributes": [
            {
                "cap_color": None,
                "class": "Model_Capper_State",
                "display_style": None,
                "surface": (0, 0),
                "version": 1,
            }
        ],
        "cap_color": None,
        "cap_offset": 0.01,
        "class": "Caps_State",
        "default_cap_offset": 0.01,
        "mesh_style": False,
        "shown": True,
        "subdivision_factor": 1.0,
        "version": 1,
    }
    import SurfaceCap.session

    SurfaceCap.session.restore_cap_attributes(cap_attributes)


registerAfterModelsCB(restore_cap_attributes)


def restore_volume_data():
    volume_data_state = {
        "class": "Volume_Manager_State",
        "data_and_regions_state": [],
        "version": 2,
    }
    from VolumeViewer import session

    session.restore_volume_data_state(volume_data_state)


try:
    restore_volume_data()
except:
    reportRestoreError("Error restoring volume data")

geomData = {"AxisManager": {}, "CentroidManager": {}, "PlaneManager": {}}

try:
    from StructMeasure.Geometry import geomManager

    geomManager._restoreSession(geomData)
except:
    reportRestoreError("Error restoring geometry objects in session")


def restoreSession_RibbonStyleEditor():
    import RibbonStyleEditor
    import SimpleSession

    userScalings = []
    userXSections = []
    userResidueClasses = []
    residueData = [
        (1, "Chimera default", "rounded", u"unknown"),
        (2, "Chimera default", "rounded", u"unknown"),
    ]
    flags = RibbonStyleEditor.NucleicDefault1
    SimpleSession.registerAfterModelsCB(
        RibbonStyleEditor.restoreState,
        (userScalings, userXSections, userResidueClasses, residueData, flags),
    )


try:
    restoreSession_RibbonStyleEditor()
except:
    reportRestoreError("Error restoring RibbonStyleEditor state")

trPickle = "gAJjQW5pbWF0ZS5UcmFuc2l0aW9ucwpUcmFuc2l0aW9ucwpxASmBcQJ9cQMoVQxjdXN0b21fc2NlbmVxBGNBbmltYXRlLlRyYW5zaXRpb24KVHJhbnNpdGlvbgpxBSmBcQZ9cQcoVQZmcmFtZXNxCEsBVQ1kaXNjcmV0ZUZyYW1lcQlLAVUKcHJvcGVydGllc3EKXXELVQNhbGxxDGFVBG5hbWVxDWgEVQRtb2RlcQ5VBmxpbmVhcnEPdWJVCGtleWZyYW1lcRBoBSmBcRF9cRIoaAhLFGgJSwFoCl1xE2gMYWgNaBBoDmgPdWJVBXNjZW5lcRRoBSmBcRV9cRYoaAhLAWgJSwFoCl1xF2gMYWgNaBRoDmgPdWJ1Yi4="
scPickle = "gAJjQW5pbWF0ZS5TY2VuZXMKU2NlbmVzCnEBKYFxAn1xA1UHbWFwX2lkc3EEfXNiLg=="
kfPickle = (
    "gAJjQW5pbWF0ZS5LZXlmcmFtZXMKS2V5ZnJhbWVzCnEBKYFxAn1xA1UHZW50cmllc3EEXXEFc2Iu"
)


def restoreAnimation():
    "A method to unpickle and restore animation objects"
    # Scenes must be unpickled after restoring transitions, because each
    # scene links to a 'scene' transition. Likewise, keyframes must be
    # unpickled after restoring scenes, because each keyframe links to a scene.
    # The unpickle process is left to the restore* functions, it's
    # important that it doesn't happen prior to calling those functions.
    import SimpleSession
    from Animate.Session import restoreKeyframes, restoreScenes, restoreTransitions

    SimpleSession.registerAfterModelsCB(restoreTransitions, trPickle)
    SimpleSession.registerAfterModelsCB(restoreScenes, scPickle)
    SimpleSession.registerAfterModelsCB(restoreKeyframes, kfPickle)


try:
    restoreAnimation()
except:
    reportRestoreError("Error in Animate.Session")


def restoreLightController():
    import Lighting

    Lighting._setFromParams(
        {
            "ratio": 1.25,
            "brightness": 1.16,
            "material": [30.0, (0.85, 0.85, 0.85), 0.0],
            "back": [
                (0.3574067443365933, 0.6604015517481455, -0.6604015517481456),
                (1.0, 1.0, 1.0),
                0.0,
            ],
            "mode": "two-point",
            "key": [
                (-0.3574067443365933, 0.6604015517481455, 0.6604015517481456),
                (1.0, 1.0, 1.0),
                1.0,
            ],
            "contrast": 0.83,
            "fill": [
                (0.2505628070857316, 0.2505628070857316, 0.9351131265310294),
                (1.0, 1.0, 1.0),
                0.0,
            ],
        }
    )


try:
    restoreLightController()
except:
    reportRestoreError("Error restoring lighting parameters")


def restoreRemainder():
    from SimpleSession.versions.v65 import (
        restoreFontInfo,
        restoreModelClip,
        restoreOpenModelsAttrs,
        restoreOpenStates,
        restoreSelections,
        restoreSilhouettes,
        restoreWindowSize,
    )

    curSelIds = []
    savedSels = []
    openModelsAttrs = {"cofrMethod": 4}
    windowSize = (888, 678)
    xformMap = {
        0: (
            ((-0.11624158307307, -0.87318323267119, 0.47332751509558), 18.779158484753),
            (0.50892550362133, -0.03438019354708, 0.061560056363412),
            True,
        )
    }
    fontInfo = {"face": ("Sans Serif", "Normal", 16)}
    clipPlaneInfo = {}
    silhouettes = {0: True, 264: True, 265: True}

    replyobj.status("Restoring window...", blankAfter=0, secondary=True)
    restoreWindowSize(windowSize)
    replyobj.status("Restoring open states...", blankAfter=0, secondary=True)
    restoreOpenStates(xformMap)
    replyobj.status("Restoring font info...", blankAfter=0, secondary=True)
    restoreFontInfo(fontInfo)
    replyobj.status("Restoring selections...", blankAfter=0, secondary=True)
    restoreSelections(curSelIds, savedSels)
    replyobj.status("Restoring openModel attributes...", blankAfter=0, secondary=True)
    restoreOpenModelsAttrs(openModelsAttrs)
    replyobj.status("Restoring model clipping...", blankAfter=0, secondary=True)
    restoreModelClip(clipPlaneInfo)
    replyobj.status("Restoring per-model silhouettes...", blankAfter=0, secondary=True)
    restoreSilhouettes(silhouettes)

    replyobj.status(
        "Restoring remaining extension info...", blankAfter=0, secondary=True
    )


try:
    restoreRemainder()
except:
    reportRestoreError("Error restoring post-model state")
from SimpleSession.versions.v65 import makeAfterModelsCBs

makeAfterModelsCBs()

from SimpleSession.versions.v65 import endRestore

replyobj.status("Finishing restore...", blankAfter=0, secondary=True)
endRestore({})
replyobj.status("", secondary=True)
replyobj.status("Restore finished.")
