{ pkgs ? import <nixpkgs> {} }:
pkgs.mkShell {
packages = [
    pkgs.bazel_5
    pkgs.jdk11
    pkgs.clang
];

inputsFrom = [ pkgs.hello pkgs.gnutar ];
}